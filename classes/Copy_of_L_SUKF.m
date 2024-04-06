classdef L_SUKF < handle
    % Unscented Kalman Filter implementation based on the L-SUKF variant
    % from the article "An In-Flight Alignment Method for GPS-assisted ..."
%     properties (Access = private)
    properties (Access = public)
        TS = 0.01;            % Sampling interval, default - 10 msec. (from xSens 100 Hz logs)
        L = 12;               % Number of the state vector parameters
        
        % For the main article
        % J = L;                % Input dimension for Sigma points
        % I = J + 2;            % Number of Sigma points
        % X_ij_Sigma_points = zeros(J, I); % ij Sigma point values
        % X_i_Sigma_points = zeros(I, 1); % i Sigma point values
        % W_m_i = zeros(L,1);   % The state weight coefficient values
        % W_c_i = zeros(L,1);   % The weight coefficients of the Sigma sampling points                       
        % W0 = 0.5;             % The initial state weight coefficient value, 0 <= W_0 < 1
        % W1 = 0.5;             % The initial state weight coefficient value
        
        alpha_coef = 0.003;       % The positive zooming factor, 0.0001 <= alpha <= 1
        beta_coef = 2;            % The coefficient to incorporate the prior knowledge of the X parameters pdfs, beta = 2 for Gaussian PDF
        k_coef = 0;               % The secondary scaling parameter
        lambda;                   % Scaling parameter
        
        W_m_i;                    % The state weight coefficient values
        W_c_i;                    % The weight coefficients of the Sigma sampling points
        
        w_b_a;                    % The Accelerometer measurement noise variance
        w_b_g;                    % The Gyroscope measurement noise variance
        w_g_n;                    % The GPS measurement noise variance
    end
    
    properties (Access = public)
        phi_b_k   = zeros(3,1);   % The tracking error of the attitude vector
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
            if size(phi_b_0) == [1,3], obj.phi_b_k = phi_b_0'; 
            else, obj.phi_b_k = phi_b_0; end
            
            if size(dV_ib_0) == [1,3], obj.dV_ib_k = dV_ib_0'; 
            else, obj.dV_ib_k = dV_ib_0; end
            
            if size(eps_b_0) == [1,3], obj.eps_b_k = eps_b_0'; 
            else, obj.eps_b_k = eps_b_0; end
            
            if size(delta_b_0) == [1,3], obj.delta_b_k = delta_b_0'; 
            else, obj.delta_b_k = delta_b_0; end
            
            % Establish value of the Scaling parameter
            obj.lambda   = (obj.alpha_coef^2) * (obj.L + obj.k_coef) - obj.L;
            % Initialize the state weight coefficient values
            obj.W_m_i    = zeros(2*obj.L + 1, 1);
            % Initialize the weight coefficients of the Sigma sampling points
            obj.W_c_i    = zeros(2*obj.L + 1, 1);
           
            obj.CalcW_m_i_coefs_way2();
            obj.CalcW_c_i_coefs_way2();
            
            obj.B        = [obj.dV_ib_k; obj.eps_b_k; obj.delta_b_k];
            obj.X_k_1    = [obj.phi_b_k; obj.B];
            obj.X_k_pred = zeros(size(obj.X_k_1));
            obj.X_k      = zeros(size(obj.X_k_1));
            
            obj.P_k_1    = P0;
            obj.P_k_pred = zeros(size(obj.P_k_1));
            obj.P_k      = zeros(size(obj.P_k_1));
            
            obj.F        = ones(size(obj.X_k, 1));            
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
%             phi_b_k_upd  = exp(dX_k(1:size(obj.phi_b_k,1), 1)) .* ...
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
            obj.G(1:3, 4:6) = R_ib_b;
            obj.G(4:6, 1:3) = -R_ib_b;            
            output = obj.G * [ones(3,1)*obj.w_b_g; ones(3,1)*obj.w_b_a; zeros(6,1)];
        end
        
        
        %% Calculate matrix R
        function output = calcR(obj)
            output = eye(3) * obj.w_g_n;
        end
        
        
        %% Calculate P(k|k-1) matrix from X(k|k-1) vector and Y Sigma points
        function output = calc_P_pred(obj, X_k_values, sigmaYpoints)
            % Calculate nu_i first (size L x (2*L+1))
            nu_i = sigmaYpoints - X_k_values;
            % Calculate the P(k|k-1)
            tempP = zeros(size(obj.P_k));
            for i=1:size(sigmaYpoints,2)
                tempP = tempP + obj.W_c_i(i,1) * nu_i(:,i) * nu_i(:,i)';
            end
            % Add the model Noise matrix
            output = tempP + obj.Q;
        end
        
        
        %% Restore the X(k|k-1) vector from Y Sigma points
        function output = calc_X_pred_from_SigmaY_points(obj, sigmaYpoints)
            phiLen  = size(obj.phi_b_k, 1);
            % Get phi^b(k|k-1)
            ksiPhi  = sigmaYpoints(1:phiLen, 1:size(sigmaYpoints,2));
            tempPhi = zeros(size(obj.phi_b_k));
            for i=1:size(sigmaYpoints,2)
                tempPhi = tempPhi + obj.W_m_i(i,1) * ksiPhi(:,i);
            end
            
            % Get B(k|k-1)
            ksiB   = sigmaYpoints(phiLen+1:phiLen+size(obj.B,1), 1:size(sigmaYpoints,2));
            tempB   = zeros(size(obj.B));
            for i=1:size(sigmaYpoints,2)
                tempB = tempB + obj.W_m_i(i,1)*ksiB(:,i);
            end
            
            output = [tempPhi; tempB];
        end
                

        %% Calculate SigmaY values as SigmaY = f(SigmaX)
        function output = calc_SigmaY_values(obj, acc, R_ib_b)
            % 1a - Calculate SigmaX points based on the X(k-1) state (size = [L, 2L+1])
            sigmaX_points = obj.Calc_Sigma_points_way2(obj.X_k_1, obj.P_k_1);
            
            % Calc A matrix values
            acc_abs = sqrt(sum(acc.^2));
            p_vect = acc / acc_abs;            
            exp_f_sks_dt = cos(obj.TS * acc_abs * ones(3)) + ...
                           (1 - cos(obj.TS * acc_abs)) * (p_vect * p_vect') + ...
                           sin(obj.TS * acc_abs * obj.getSkewSym(p_vect));
            exp_R_ib_b_dt = exp(obj.TS * R_ib_b);
            exp_minus_R_ib_b_dt = exp(-obj.TS * R_ib_b);
            
            obj.F(1:3, 4:6)   = exp_f_sks_dt;
            obj.F(1:3, 10:12) = exp_R_ib_b_dt;
            obj.F(4:6, 7:9)   = exp_minus_R_ib_b_dt;
            
            output = obj.F * sigmaX_points;
            
%             obj.G(1:3, 4:6)   = R_ib_b;
%             obj.G(4:6, 1:3)   = -R_ib_b;
%             output = obj.F * (sigmaX_points + obj.G * [obj.w_b_g; obj.w_b_a]);            
        end
            
        
        %% Following the article "UKF for non-linear estimation"
        function output = Calc_Sigma_points_way2(obj, mean_X, Px_matrix)
            X_i_Sigma_points = zeros(obj.L, 2*obj.L+1);
            
            % Calc the second term of the SigmaX point formula
            P_sq_root = sqrt((obj.L + obj.lambda) * Px_matrix);
            % i = 0
            X_i_Sigma_points(:,1) = mean_X;
            % for i = 1..L
            for i = 2 : obj.L+1
                X_i_Sigma_points(:,i) = mean_X + P_sq_root(i-1,:)';
            end
            % for i = L+1...2L
            for i = obj.L+2 : 2*obj.L+1
                X_i_Sigma_points(:,i) = mean_X - P_sq_root(i-(obj.L+1),:)';
            end
            output = X_i_Sigma_points;
        end
        
        
        %% Calculate the W_m_i coefficients
        function CalcW_m_i_coefs_way2(obj)
            % W_m_0
            obj.W_m_i(1,1) = obj.lambda / (obj.L + obj.lambda);
            % Calculate all remaining values of W_m_i, i=1..2L
            for i = 2 : 2*obj.L+1
                obj.W_m_i(i,1) = 1 / 2 / (obj.L + obj.lambda);
            end
        end
        
        
        %% Calculate the W_m_i coefficients
        function CalcW_c_i_coefs_way2(obj)
            % W_c_0
            obj.W_c_i(1,1) = obj.lambda / ((obj.L + obj.lambda) + ...
                (1 - obj.alpha_coef^2 + obj.beta_coef));
            % Calculate all remaining values of W_c_i, i=1..2L
            for i = 2 : 2*obj.L+1
                obj.W_c_i(i,1) = 1 / 2 / (obj.L + obj.lambda);
            end
        end
        
        
        %% Following the initial_article
%         function Calc_Sigma_points_way1(obj, mean_X, P_matrix)
%             % Fill in for i=0
%             curr_ksi = zeros(obj.J, 1);
%             obj.X_ij_Sigma_points(:,1) = curr_ksi;
%             % Fill in for i=1,3,5,...
%             curr_ksi = ones(obj.J, 1) * (-1 / sqrt(2 * obj.W_m_i(2,1)));
%             obj.X_ij_Sigma_points(:,2:2:obj.I-1) = curr_ksi;
%             % Fill in for i=2,4,6,...
%             curr_ksi = ones(obj.J, 1) * (-1 / sqrt(2 * obj.W_m_i(2,1)));
%             curr_ksi(1,1) = 1 / sqrt(2 * obj.W_m_i(2,1));
%             obj.X_ij_Sigma_points(:,3:2:obj.I-1) = curr_ksi;
%             % Fill in for i=obj.I, i.e. J+1 if 'i' is counted from 0
%             curr_ksi = zeros(obj.J, 1);
%             curr_ksi(obj.J,1) = 1 / sqrt(2 * obj.W_m_i(2,1));
%             obj.X_ij_Sigma_points(:,obj.I) = curr_ksi;
%             
%             % Get the sqrt(P_matrix); it should be J x J
%             sqrtP = sqrt(P_matrix);
%             
%             % Calculate the 'i' Sigma points
%             for i=1:obj.I
%                 obj.X_i_Sigma_points(i,1) = mean_X(i,1) + ...   % ERRORR!!!!
%                     obj.alpha * sqrtP(i,:) * obj.X_ij_Sigma_points(:,i);
%             end
%         end
        
        %% Calculate the W_m_i coefficients
%         function CalcW_m_i_coefs_way1(obj)
%             % W_m_0
%             obj.W_m_i(1,1) = obj.W0^2 / obj.alpha_coef^2 + (1 - 1/obj.alpha_coef^2);
%             % W_m_1
%             obj.W_m_i(2,1) = (1 - obj.W0) / (2^(obj.L-1) * obj.alpha_coef^2);
%             % W_m_2
%             obj.W_m_i(3,1) = obj.W_m_i(2,1);
%             % Calculate all remaining values of W_m_i, i=3..n+1
%             for i=4:obj.L+1
%                 obj.W_m_i(i,1) = 2^(i-1-2) * obj.W1 / obj.alpha_coef^2;
%             end
%         end
        
        %% Calculate the W_m_i coefficients
%         function CalcW_c_i_coefs_way1(obj)
%             % W_c_0
%             obj.W_c_i(1,1) = obj.W_m_i(1,1) + (obj.W0 + 1 + obj.beta_coef - obj.alpha_coef^2);
%             % Calculate all remaining values of W_c_i, i=2..n+1
%             for i=4:obj.L+1
%                 obj.W_c_i(i,1) = obj.W_m_i(i,1);
%             end
%         end
        
        
        %% Method provides Skew-Symmetric matrix from the input vector
        function output = getSkewSym(obj, inpVect)
            output = zeros(3);
            if length(inpVect) == 3
                output = [      0     -inpVect(3)  inpVect(2);...
                           inpVect(3)      0      -inpVect(1);...
                           -inpVect(2) inpVect(1)      0         ];
            else
                warning('L_SUKF_Method.getSkewSym(): Wrong input vector dimension!');
            end
        end

    end
end








