classdef L_SUKF < handle
    % Unscented Kalman Filter implementation based on the L-SUKF variant
    % from the article "An In-Flight Alignment Method for GPS-assisted ..."
%     properties (Access = private)
    properties (Access = public)
        TS = 0.01;                % Sampling interval, default - 10 msec. (from xSens 100 Hz logs)
        L  = 12;                  % Number of the state vector parameters
        termsN = 3;               % Number of terms for series expansion for SO(3) to so(3)
        
        X_ij_Sigma_points;        % ij Sigma point values
        X_i_Sigma_points;         % i Sigma point values                   
        W0 = 0.5;                 % The initial state weight coefficient value, 0 <= W_0 < 1
        W_i;                      % The weight coefficient values, W_i, 0 <= i <= L+1
        
        alpha_coef = 0.003;       % The positive zooming factor, 0.0001 <= alpha <= 1
        beta_coef = 2;            % The coefficient to incorporate the prior knowledge of the X parameters pdfs, beta = 2 for Gaussian PDF
        
        W_m_i;                    % The state weight coefficient values
        W_c_i;                    % The weight coefficients of the Sigma sampling points
        
        w_b_a;                    % The Accelerometer measurement noise variance
        w_b_g;                    % The Gyroscope measurement noise variance
        w_g_n;                    % The GPS measurement noise variance
    end
    
    properties (Access = public)
        phi_b_k   = zeros(3,1);   % The tracking error of the attitude vector
        phi_b_k_SO3 = zeros(3);   % The phi_b_k vector in the SO(3) space
        dV_ib_k   = zeros(3,1);   % The observation vector error        
        eps_b_k   = zeros(3,1);   % The Gyroscope bias
        delta_b_k = zeros(3,1);   % The Accelerometer bias
        
        % State vector X = [phi_b_k; B], where B = [dV_ib_k; eps_b_k; delta_b_k]        
        X_k;                      % UKF State vector at the current step (k)
        X_k_pred;                 % Predicted UKF State vector for the next step k
        X_k_1;                    % UKF State vector at the previous step (k-1)
        B;                        % Matrix B = [dV_ib_k; eps_b_k; delta_b_k], part of matrix X
        
        P_k;                      % P covariance matrix at the current step (k)
        P_k_pred;                 % Predicted P covariance matrix for the next step k
        P_k_1;                    % P covariance matrix at the previous step (k-1)
        Q;                        % Process noise covariance matrix (sizeof(X_k) x sizeof(X_k))
        
        F;                        % The matrix F, as matrix exp(dt * matrix_A)
        G;                        % The matrix G, as Q = G*W       
        K;                        % Kalman Gain matrix 
        H;                        % H matrix
        Z;                        % Measurement vector
        Z_pred;                   % Predicted Measurement vector
        R;                        % Measurement noise covariance matrix (sizeof(Z) x sizeof(Z))
    end
    
    
    
    methods  (Access = public)        
        %% Constructor
        function obj = L_SUKF(varargin)
            for i = 1:2:nargin
                if     strcmpi(varargin{i}, 'TS'), obj.setTS(varargin{i+1});
                elseif strcmpi(varargin{i}, 'alpha-coef'), obj.setAlphaCoef(varargin{i+1});
                elseif strcmpi(varargin{i}, 'beta-coef'), obj.setBetaCoef(varargin{i+1});
                elseif strcmpi(varargin{i}, 'acc-noise'), obj.set_w_b_a(varargin{i+1});
                elseif strcmpi(varargin{i}, 'gyro-noise'), obj.set_w_b_g(varargin{i+1});
                elseif strcmpi(varargin{i}, 'gps-noise'), obj.set_w_g_n(varargin{i+1});
                else
                    error('L_SUKF.m: Invalid argument!');
                end
            end            
            % Necessary to define here because it is used later in Init-method
            obj.X_k = zeros(obj.L,1);
        end
        
        
        %% Function for initialization of the first sample of the method output
        function obj = Init(obj, phi_b_0, dV_ib_0, eps_b_0, delta_b_0, P0)
            if size(phi_b_0,1) == 1, obj.phi_b_k = phi_b_0'; 
            else, obj.phi_b_k = phi_b_0; end
            
            if size(dV_ib_0,1) == 1, obj.dV_ib_k = dV_ib_0'; 
            else, obj.dV_ib_k = dV_ib_0; end
            
            if size(eps_b_0,1) == 1, obj.eps_b_k = eps_b_0'; 
            else, obj.eps_b_k = eps_b_0; end
            
            if size(delta_b_0,1) == 1, obj.delta_b_k = delta_b_0'; 
            else, obj.delta_b_k = delta_b_0; end
            
            % Initialize the state weight coefficient values
            obj.W_m_i    = zeros(obj.L+2, 1);
            % Initialize the weight coefficients of the Sigma sampling points
            obj.W_c_i    = zeros(obj.L+2, 1);
            
            obj.X_ij_Sigma_points = zeros(obj.L, obj.L+2);
            obj.X_i_Sigma_points  = zeros(obj.L, obj.L+2);
            
            %obj.Calc_Xij_coefficient_matrix_incor();
            obj.X_ij_Sigma_points = obj.Calc_Xij_coefficient_matrix_ver2();
            obj.CalcW_m_i_coefs();
            obj.CalcW_c_i_coefs();
            
            obj.B        = [obj.dV_ib_k; obj.eps_b_k; obj.delta_b_k];
            obj.X_k_1    = [obj.phi_b_k; obj.B];
            obj.X_k_pred = zeros(size(obj.X_k_1));
            obj.X_k      = zeros(size(obj.X_k_1));
            
            obj.P_k_1    = P0;
            obj.P_k_pred = zeros(size(obj.P_k_1));
            obj.P_k      = zeros(size(obj.P_k_1));
            
            obj.F        = zeros(size(obj.X_k,1));            
            obj.G        = zeros(size(obj.X_k,1));
            
            dV_size      = size(obj.dV_ib_k,1);
            obj.H        = [eye(dV_size) zeros(dV_size, size(obj.X_k,1)-dV_size)];
            obj.Z        = zeros(size(dV_size,1));
            
            obj.R        = obj.calcR();
        end
        
        
        %% Main function where sensor measurements (Accel) are passed to the filter
        function obj = Update(obj, acc, R_ib_b, V_ib_k, V_in_k, R_ib_in)
            % Initialize necessary vars
            obj.X_k_pred = zeros(size(obj.X_k_1));
            obj.X_k      = zeros(size(obj.X_k_1));
            obj.P_k_pred = zeros(size(obj.P_k_1));
            obj.P_k      = zeros(size(obj.P_k_1));
            
            % Stage 1: Prediction - Update KF State Update Equation
            %       1a:Predict Y Sigma point values (for the X(k) state)
            sigmaYpoints = obj.calc_SigmaY_values(acc, R_ib_b);
            %       1b:Restore the X(k|k-1) vector from Y Sigma points
            obj.X_k_pred = obj.calc_X_pred_from_SigmaY_points(sigmaYpoints);
            %       1c:Calculate P(k|k-1) matrix from X(k|k-1) vector and Y Sigma points
            obj.Q        = obj.calcQ(R_ib_b);
            %       1d:Calculate P(k|k-1) matrix from X(k|k-1) vector and Y Sigma points
            obj.P_k_pred = obj.calc_P_pred(obj.X_k_pred, sigmaYpoints);
            
            
            % Stage 2: Correction - Update KF State Update Equation
            %       2a:Calculate measurement and predicted measurement matrices
            obj.Z        = V_ib_k - R_ib_in * V_in_k;
            obj.Z_pred   = obj.H * obj.X_k_pred;
            
            %       2b:Calculate Pzz(k) and Pxz(k)
            P_zk_zk      = obj.H * obj.P_k_pred * obj.H' + obj.R;
            P_xk_zk      = obj.P_k_pred * obj.H';
            
            %       2c:Calculate Filter gain matrix K(k)
            obj.K        = P_xk_zk / P_zk_zk;
            
            %       2d:Calculate the difference of X vector values
            dX_k         = obj.K * (obj.Z - obj.Z_pred);
            
            %       2e:Update the both parts of the X vector
%             phi_b_k_upd  = exp(dX_k(1:size(obj.phi_b_k_SO3,1), 1)) .* ...
%                               obj.X_k_pred(1:size(obj.phi_b_k,1), 1);
            phi_b_k_upd  = dX_k(1:size(obj.phi_b_k,1), 1) + ...
                              obj.X_k_pred(1:size(obj.phi_b_k,1), 1);
            B_upd        = dX_k(size(obj.phi_b_k,1)+1:end, 1) + ...
                              obj.X_k_pred(size(obj.phi_b_k,1)+1:end, 1);
            
            obj.phi_b_k  = phi_b_k_upd;
            obj.B        = B_upd;
            
            % Note that B = [dV_ib_k; eps_b_k; delta_b_k]
            obj.dV_ib_k  = B_upd(1:size(obj.dV_ib_k,1), 1);
            obj.eps_b_k  = B_upd(size(obj.dV_ib_k,1)+1: ...
                                 size(obj.dV_ib_k,1)+size(obj.eps_b_k,1), 1);
            obj.delta_b_k= B_upd(size(obj.dV_ib_k,1)+size(obj.eps_b_k,1)+1:...
                                 size(obj.dV_ib_k,1)+size(obj.eps_b_k,1)+size(obj.delta_b_k,1), 1);
            obj.X_k      = [obj.phi_b_k; obj.B];
            
            %       2f:Update the Covariance matrix (P(k)) of the system
            obj.P_k      = obj.P_k_pred - obj.K * P_zk_zk * obj.K';
            
            
            % Update the X(k-1) matrix for the next iteration
            obj.X_k_1    = obj.X_k;
            obj.P_k_1    = obj.P_k;
        end
        
        
        
        %% Setter for sampling interval value
        function setTS(obj, newTS)
            obj.TS = newTS;
        end
        %% Getter for sampling interval value
        function output = getTS(obj)
            output = obj.TS;
        end
        
        %% Setter for Process noise covariance matrix
        function setQ(obj, newSQ)
            obj.Q = newSQ;
        end
        %% Getter for Process noise covariance matrix
        function output = getQ(obj)
            output = obj.Q;
        end
        
        %% Setter for Measurement noise covariance matrix
        function setR(obj, newSR)
            obj.R = newSR;
        end
        %% Getter for Measurement noise covariance matrix
        function output = getR(obj)
            output = obj.R;
        end
        
        %% Setter for the positive zooming factor
        function setAlphaCoef(obj, newAlphaCoef)
            obj.alpha_coef = newAlphaCoef;
        end
        %% Getter for the positive zooming factor
        function output = getAlphaCoef(obj)
            output = obj.alpha_coef;
        end
        
        %% Setter for the Beta coefficient
        function setBetaCoef(obj, newBetaCoef)
            obj.beta_coef = newBetaCoef;
        end
        %% Getter for the Beta coefficient
        function output = getBetaCoef(obj)
            output = obj.beta_coef;
        end
   
        
        %% Setter for the Accel measument noise variance value
        function set_w_b_a(obj, new_w_b_a)
            obj.w_b_a = new_w_b_a;
        end
        %% Getter for the Accel measument noise variance value
        function output = get_w_b_a(obj)
            output = obj.w_b_a;
        end
        
        %% Setter for the Gyro measument noise variance value
        function set_w_b_g(obj, new_w_b_g)
            obj.w_b_g = new_w_b_g;
        end
        %% Getter for the Gyro measument noise variance value
        function output = get_w_b_g(obj)
            output = obj.w_b_g;
        end
        
        %% Setter for the GPS measument noise variance value
        function set_w_g_n(obj, new_w_g_n)
            obj.w_g_n = new_w_g_n;
        end
        %% Getter for the GPS measument noise variance value
        function output = get_w_g_n(obj)
            output = obj.w_g_n;
        end
    end    
    
    
    methods (Access = private) 
        
        %% Calculate matrix Q
        function output = calcQ(obj, R_ib_b)
            obj.G(1:3, 1:3) = -R_ib_b;
            obj.G(4:6, 4:6) = R_ib_b;            
            output = obj.TS * obj.G * [ones(3,1)*obj.w_b_g; ones(3,1)*obj.w_b_a; zeros(6,1)];
        end
        
        
        %% Calculate matrix R
        function output = calcR(obj)
            output = eye(3) * obj.w_g_n;
        end
        
        
        %% Calculate P(k|k-1) matrix from X(k|k-1) vector and Y Sigma points
        function output = calc_P_pred(obj, X_k_values, sigmaYpoints)
            phiLen  = size(obj.phi_b_k, 1);
            ksiB_sigma_points = sigmaYpoints(phiLen+1:phiLen+size(obj.B,1), :);
            
            % Calculate nu_a
            inv_phi_b_k = inv(obj.phi_b_k_SO3);
            nu_a = zeros(phiLen, size(sigmaYpoints,2));
            for i=1:size(sigmaYpoints,2)
                temp = obj.so3_to_SO3(obj.getSkewSym(sigmaYpoints(1:phiLen, i))) * inv_phi_b_k;
                nu_a(:,i) = obj.so3_to_euclid(obj.SO3_to_so3(temp));
            end
            
            % Calculate nu_b
            nu_b = zeros(size(obj.B,1), size(sigmaYpoints,2));
            for i=1:size(sigmaYpoints,2)
                nu_b(:,i) = ksiB_sigma_points(:,i) - X_k_values(phiLen+1:phiLen+size(obj.B,1), 1);
            end
            
            % Combine nu-matrices
            nu_comb = [nu_a; nu_b];
            
            % Calculate the P(k|k-1)
            tempP = zeros(size(obj.P_k));
            for i=1:size(sigmaYpoints,2)
                tempP = tempP + obj.W_c_i(i,1) * nu_comb(:,i) * (nu_comb(:,i)');
            end
            
            % Add the model Noise matrix
            output = tempP + obj.Q;
        end
        
        
        %% Restore the X(k|k-1) vector from Y Sigma points
        function output = calc_X_pred_from_SigmaY_points(obj, sigmaYpoints)
            phiLen  = size(obj.phi_b_k, 1);
            % Get phi^b(k|k-1) in SO(3) space for Pk calculation
            obj.phi_b_k_SO3 = obj.Calc_phi_b_pred_SO3_Alg1(sigmaYpoints);
            
            % Get phi^b(k|k-1) (not written in the paper)
            ksiPhi = obj.so3_to_euclid(obj.SO3_to_so3(obj.phi_b_k_SO3));
            
            % Get B(k|k-1)
            ksiB    = sigmaYpoints(phiLen+1:phiLen+size(obj.B,1), :);
            tempB   = zeros(size(obj.B));
            for i=1:size(sigmaYpoints,2)
                tempB = tempB + obj.W_m_i(i,1)*ksiB(:,i);
            end
            
            output = [ksiPhi; tempB];
        end
                

        %% Calculate SigmaY values as SigmaY = f(SigmaX)
        function output = calc_SigmaY_values(obj, acc, R_ib_b)
            % Calculate SigmaX points based on the X(k-1) state (size = [L, 2L+1])
            sigmaX_points = obj.Calc_Ksi_i_Sigma_points(obj.X_k_1, obj.P_k_1);
            
            % Calculate matrix F
            obj.F(1:3, 7:9)   = -R_ib_b;
            obj.F(4:6, 1:3)   = obj.getSkewSym(acc);
            obj.F(4:6, 10:12) = R_ib_b;
            
            % Compose State Prediction Equation as X(k|k-1) = (I+dT*F)* X(k-1)
            output = (eye(obj.L) + obj.TS * obj.F) * sigmaX_points;
        end
            
        
        %% Following the initial_article
        function output = Calc_Ksi_i_Sigma_points(obj, mean_X, P_matrix)
            output = zeros(obj.L, obj.L+2);
            
            % Calc the second term of the SigmaX point formula
            P_sq_root = chol(P_matrix);
            %P_sq_root = sqrtm(P_matrix);
            
            % The loop for main calculation of the Sigma X points
            for i=1:obj.L+2
                output(:,i) = mean_X + obj.alpha_coef * P_sq_root * obj.X_ij_Sigma_points(:,i);
            end
        end


        %% Fill in the Xij-coefficient matrix
        %    Min skew simplex UT approach (from the article in the file - 
        %      "Scale_corrected_minimal_skew_simplex_sampling_UKF_for_BLDCM_sensorless.pdf")
        function output = Calc_Xij_coefficient_matrix_ver2(obj)
            output = nan(obj.L, obj.L+2);
            
            % Calculate the W_i coefficients for ksi^j_i formulas
            obj.W_i = zeros(obj.L+1,1);
            obj.W_i(1,1) = (1 - obj.W0)/2^obj.L;
            obj.W_i(2,1) = obj.W_i(1);
            for i=3:obj.L+1
                obj.W_i(i,1) = 2^(i-1) * obj.W_i(1,1);
            end
            
            % For i=0,1,...,L+1
            for i=0:obj.L+1
                output(:,i+1) = get_Xij_vector(obj, obj.L, i);
            end
        end

        
        function output = get_Xij_vector(obj, j, i)
            if i == 0
                if j == 1
                    output = 0;
                    return
                else
                    output = [obj.get_Xij_vector(j-1, 0); 0];
                    return
                end
            elseif i == 1 || i == 2
                if j == 1
                    if i == 1
                        output = -1 / sqrt(2 * obj.W_i(1,1));
                        return
                    end
                    if i == 2
                        output = 1 / sqrt(2 * obj.W_i(1,1));
                        return 
                    end
                else
                    output = [obj.get_Xij_vector(j-1, i); -1 / sqrt(2 * obj.W_i(j+1,1))];
                    return
                end
            elseif i > 2
                if i <= j
                    output = [obj.get_Xij_vector(j-1, i); -1 / sqrt(2 * obj.W_i(j+1,1))];
                    return
                elseif i == (j + 1)
                    output = [zeros(j-1,1); 1 / sqrt(2 * obj.W_i(j+1,1))];
                    return
                end
            end
                
        end
            
        
        %% Calculate the W_m_i coefficients
        function CalcW_m_i_coefs(obj)
            % W_m_0
            obj.W_m_i(1,1) = obj.W0/(obj.alpha_coef^2) + (1 - 1/(obj.alpha_coef^2));
            % W_m_1
            obj.W_m_i(2,1) = (1 - obj.W0) / (2^obj.L) / (obj.alpha_coef^2);
            % W_m_2
            obj.W_m_i(3,1) = obj.W_m_i(2,1);
            % Calculate all remaining values of W_m_i, i=3..n+1
            for i=4:obj.L+2
                obj.W_m_i(i,1) = 2^(i-1-2) * obj.W_i(1,1) / obj.alpha_coef^2;
            end
        end
        
        %% Calculate the W_c_i coefficients
        function CalcW_c_i_coefs(obj)
            % W_c_0
            obj.W_c_i(1,1) = obj.W_m_i(1,1) + 1 + obj.beta_coef - obj.alpha_coef^2;
            % Calculate all remaining values of W_c_i, i=2..n+1
            for i=2:obj.L+2
                obj.W_c_i(i,1) = obj.W_m_i(i,1);
            end
        end
        
        
        %% Algorithm 1 (from the main article) implementation
        function output = Calc_phi_b_pred_SO3_Alg1(obj, sigmaYPoints)
            % Step 0 - Get all Zi matrices from the (1:3,1) values of each
            %            Sigma point except the 13th one
            matr_Zi = sigmaYPoints(1:3, 1:obj.L+1);
            
            % Step 1 - H = Z0 (convert ksi0 vector to SO3 matrix)
            matr_H = obj.so3_to_SO3(obj.getSkewSym(matr_Zi(:,1)));
            
            % Step 2 - Calc Omega matrix
            omega_matr = zeros(3,1);
            for i=2:obj.L+1
                temp = obj.SO3_to_so3(obj.so3_to_SO3(obj.getSkewSym(matr_Zi(:,i))) / matr_H);
                omega_matr = omega_matr + obj.W_m_i(i,1) * obj.so3_to_euclid(temp);
            end
            % Step 3 - Calc H matrix
            output = obj.so3_to_SO3(obj.getSkewSym(omega_matr)) * matr_H;
        end
        
        
        %% Convertion from SO(3) space to so(3) space
        function output = SO3_to_so3(obj, input_SO3_matrix)
            output = zeros(3);
            if (numel(input_SO3_matrix) == 9)
                for i=1:obj.termsN+1
                    output = output + ((-1)^(i-1) / (i-1+1)) * (input_SO3_matrix - diag(ones(3,1)))^(i-1+1);
                end
            else
                warning('L_SUKF_Method.SO3_to_so3(): Wrong input matrix dimension!');
            end
        end
        
        
        %% Convertion from so(3) space to SO(3) space
        function output = so3_to_SO3(obj, input_so3_matrix)
            output = zeros(3);
            if (numel(input_so3_matrix) == 9)
                vect_phi = obj.so3_to_euclid(input_so3_matrix);
                len_phi  = sqrt(vect_phi(1)^2 + vect_phi(2)^2 + vect_phi(3)^2);
                vect_p   = vect_phi / len_phi;
                
                output = cos(len_phi) * diag(ones(3,1)) + ...
                    (1 - cos(len_phi))* vect_p * (vect_p') + ...
                    sin(len_phi) * obj.getSkewSym(vect_p);

            else
                warning('L_SUKF_Method.so3_to_SO3(): Wrong input matrix dimension!');
            end
        end
        
        
        %% Convertion from so(3) space to Euclidian space
        function output = so3_to_euclid(obj, input_so3_matrix)
            output = zeros(3,1);
            if (numel(input_so3_matrix) == 9)
                output(1,1) = input_so3_matrix(3,2);
                output(2,1) = input_so3_matrix(1,3);
                output(3,1) = input_so3_matrix(2,1);
            else
                warning('L_SUKF_Method.so3_to_euclid(): Wrong input matrix dimension!');
            end
        end
        
        
        %% Method provides Skew-Symmetric matrix from the input vector
        function output = CalcRM(obj, anglesZYX)
            output = zeros(3);
            if length(anglesZYX) == 3
                roll  = deg2rad(anglesZYX(1));
                pitch = deg2rad(anglesZYX(2));
                yaw   = deg2rad(anglesZYX(3));
                
                output(1,1) = cos(yaw) * cos(pitch);
                output(1,2) = cos(yaw) * sin(pitch) * sin(roll) - sin(yaw) * cos(roll);
                output(1,3) = cos(yaw) * sin(pitch) * cos(roll) + sin(yaw) * sin(roll);
                
                output(2,1) = sin(yaw) * cos(pitch);
                output(2,2) = sin(yaw) * sin(pitch) * sin(roll) + cos(yaw) * cos(roll);
                output(2,3) = sin(yaw) * sin(pitch) * cos(roll) - cos(yaw) * sin(roll);
                
                output(3,1) = -sin(pitch);
                output(3,2) = cos(pitch) * sin(roll);
                output(3,3) = cos(pitch) * cos(roll);
            else
                warning('L_SUKF_Method.CalcRM(): Wrong input vector dimension!');
            end
        end
        
        
        %% Method provides Skew-Symmetric matrix from the input vector
        function output = getSkewSym(obj, inpVect)
            output = zeros(3);
            if (length(inpVect) == 3) && (numel(inpVect) == 3)
                output = [      0     -inpVect(3)  inpVect(2);...
                           inpVect(3)      0      -inpVect(1);...
                           -inpVect(2) inpVect(1)      0         ];
            else
                warning('L_SUKF_Method.getSkewSym(): Wrong input vector dimension!');
            end
        end

        
        %% Fill in the Xij-coefficient matrix
        function output = Calc_Xij_coefficient_matrix_ver1(obj)
            output = zeros(obj.L, obj.L+2);
            % For i=0, the 1st column has been already filled with zeros
            
            coefW1 = 1/(2 * obj.W_m_i(2,1))^(1/2);
            % For i=1,
            output(1,2) = (-1)^2 * coefW1;
            output(2:end,2) = -coefW1 * ones(obj.L-1,1);
            % For i=2,
            output(1,3) = (-1)^1 * coefW1;
            output(2:end,3) = -coefW1 * ones(obj.L-1,1);
            % For i=3,4,...,L+1
            for i=3:obj.L+1
                output(1,i) = (-1)^i * coefW1;
                output(2:end,i) = -coefW1 * ones(obj.L-1,1);
            end
            % For i=L+2
            output(end,end) = coefW1;
        end
        
    end
end








