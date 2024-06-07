classdef IFA_Method < handle
    % In-Flight Alignment Method based on IMU and GPS data for defining the 
    %   Rotation matrix of the high-speed rotation object
    
    properties (Access = public)        
        curr_grav = [0; 0; -9.81];% Contains the current values of the Gravity vector

        phi_b     = zeros(3,1);  % Tracking error of the RM^b(t)_ib matrix obtained from the UKF
        dV_ib_k   = zeros(3,1);  % Tracking error of the V^ib_k vector obtained from the UKF
        gBias     = zeros(3,1);  % Tracking Gyro bias values obtained from the UKF
        aBias     = zeros(3,1);  % Tracking Accel bias values obtained from the UKF
        
        RM_bt_ib  = eye(3);      % RM of the current Body position with respect to the Body Initial position
        RM_ib_in  = eye(3);      % RM of the Body initial position with respect to the Nav Frame initial position
        RM_in_nt  = eye(3);      % RM of the Initial Nav frame position with respect to the current Nav frame position
        RM_b_n    = eye(3);      % RM of the Body current position with respect to the Nav Frame current position
        
        V_ib_k    = zeros(3,1);  % Attitude determination vector in "ib" frame
        V_in_k    = zeros(3,1);  % Attitude determination vector in "in" frame
    end
    
%     properties (Access = private)
    properties (Access = public)
        TS_IMU    = 0.01;             % IMU sampling interval, default - 10 msec. (from xSens 100 Hz logs)
        TS_GPS    = 0.2;              % GPS sampling interval, default - 200 msec. (from xSens 5 Hz logs)
        AccelBias = 0.005 * 9.81;     % m/sec^2, bias of the Accel meas, default - 5mg
        AccelVar  = (9.81^2)*1e-6;    % (m/sec^2)^2, Accel meas variance, default - (1e-6 g)^2.
        GyroBias  = deg2rad(60/3600); % rad/sec, bias of the Gyro meas, default - 60 deg per hour
        GyroVar   = (deg2rad(0.017))^2; % (rad/sec)^2, Gyro meas variance, default - (0.017 deg/sec)^2
        alp_coef  = 0.97;             % Alpha-coefficient for Gravity vector extraction
        L_in      = 0;                % The carrier latitude value at the start time
        lambda_in = 0;                % The carrier longitude value at the start time
        dx_gps    = 1;                % m, GPS position error
        dv_gps    = 0.1;              % m/sec., random GPS velocity error
        
        w_ie      = 7.292115e-5;    % Const, the Earth rotation rate, rad/sec
        w_n_in    = 2e-4*ones(3,1); % Navigation frame rotation rate with respect to the initial position of this frame
        w_n_in_ss = eye(3);         % Skew-symmetric matrix for the w_n_in
        w_n_ie    = 0*ones(3,1);    % the Earth rotation rate wrt the inertial frame
        w_n_en    = 0*ones(3,1);    % the angular rate of the Navigation frame wrt the Earth frame
        
        prev_t    = 0;              % Timestamp initial value
        prev_acc  = 1e-5*ones(3,1); % Previous value (at t(k-1)) of Body Frame acceleration
        prev_gyro = 1e-5*ones(3,1); % Previous value (at t(k-1)) of Body Frame gyroscope
        prev_vel  = 1*ones(3,1);    % Previous value (at t(k-1)) of Nav Frame velocity
        
        L_SUKF_filt;                % Instance of the L-SUKF filter    
    end
    
    methods  (Access = public)        
        %% Constructor
        function obj = IFA_Method(varargin)
            for i = 1:2:nargin
                if     strcmpi(varargin{i}, 'L_in'), obj.set_L_in(varargin{i+1});
                elseif strcmpi(varargin{i}, 'lambda_in'), obj.set_lambda_in(varargin{i+1});
                elseif strcmpi(varargin{i}, 'TS_IMU'), obj.setTS_IMU(varargin{i+1});
                elseif strcmpi(varargin{i}, 'TS_GPS'), obj.setTS_GPS(varargin{i+1});
                else
                    error('IFA_Method.m: Invalid argument!');
                end
            end            
        end
        
        
        %% Function for initialization of the first sample of the method output
        function obj = Initialize(obj, init_t, initAcc, initGyro, initVel)
            % Initialize initial timestamp
            obj.prev_t = init_t;
            % Initialize the object initial Body acceleration
            obj.prev_acc = initAcc;
            % Initialize the object initial Body gyroscope
            obj.prev_gyro = initGyro;
            % Initialize the object initial velocity
            obj.prev_vel = initVel;
            % Get the skew-sym.matrix of the w_n_in
            obj.w_n_in_ss = obj.getSkewSym( obj.w_n_in );
            % Initialize the object initial Gyro bias values
            obj.gBias    = obj.GyroBias;
            % Initialize the object initial Accel bias values
            obj.aBias    = obj.AccelBias;
            
            % Create the instance of L-SUKF fitler
            obj.L_SUKF_filt = L_SUKF('TS', obj.TS_IMU, ...
                                     'acc-noise', obj.AccelVar, ...
                                     'gyro-noise', obj.GyroVar, ...
                                     'gps-noise', (obj.dx_gps + obj.dv_gps * obj.TS_IMU)^2 );
                 
            % Initialize L-SUKF fitler: init(phi_b_0, dV_ib_0, eps_b_0, delta_b_0, P0)
            obj.L_SUKF_filt.Init(zeros(3,1), ...
                                 zeros(3,1), ...
                                 obj.GyroBias * ones(3,1), ...
                                 obj.AccelBias * ones(3,1), ...
                                 diag(1 * ones(size(obj.L_SUKF_filt.X_k, 1),1)) );
        end
        
        
        %% Main method where new IMU measurements (A and G) are passed to the method
        function obj = updateIMUData(obj, curr_t, acc, gyro, newVel)
            % 1. Update Gyro part
            modGyro = gyro - obj.gBias;
            % 1. Update Accel part
            modAccel = acc - obj.aBias;
            
            obj.updateGravityVector(modAccel);
            
            %% ========== Gyro part ========== %%
            %  2. Find the RM^b(t)_ib from the differential equation
            obj.RM_bt_ib = obj.getCurrRM_bt(modGyro);
            %  3. Update the RM using the phi_b error term obtained from UKF
            obj.RM_bt_ib = (eye(3) - obj.getSkewSym(obj.phi_b)) * obj.RM_bt_ib;           
            
            %% ========== Accel part ========== %%            
            % 2. Calculate Attitude determination vector in "ib" frame
            obj.V_ib_k = obj.getV_ib_k(curr_t, modAccel, modGyro);
            % 3. Update the Attitude determination vector by the error term obtaine from the UKF
            obj.V_ib_k = obj.V_ib_k + obj.dV_ib_k;
            
            % Calculate Attitude determination vector in "in" frame
            obj.V_in_k = obj.getV_in_k(curr_t, newVel);
                        
            
            %% ========== SVD application for R_ib_in calculation ========== %%
            [U, ~, V] = svd(obj.V_ib_k * obj.V_in_k');
            obj.RM_ib_in = U * V';
            if det(obj.RM_ib_in) < 0
                V(:, end) = -V(:, end);
                obj.RM_ib_in = U * V';
                if det(obj.RM_ib_in) < 0
                    warning('IFA_Method.updateIMUData(): RM_ib_in has negative determinant!');
                end
            end
                
            %% ========== Calc current RM of Body to Nav frame ========== %%
            obj.RM_b_n = obj.RM_bt_ib * obj.RM_ib_in * obj.RM_in_nt;
            
            
            % Update the L-SUKF: update(acc, R_ib_b, V_ib_k, V_in_k, R_ib_in)
            obj.L_SUKF_filt.Update(acc, obj.RM_bt_ib', obj.V_ib_k, obj.V_in_k, obj.RM_ib_in);
            
            % Get the results from the L-SUFK and update the parameters of the main method
            obj.phi_b   = obj.L_SUKF_filt.phi_b_k;
            obj.dV_ib_k = obj.L_SUKF_filt.dV_ib_k;
            obj.gBias   = obj.L_SUKF_filt.eps_b_k;
            obj.aBias   = obj.L_SUKF_filt.delta_b_k;
            
            fprintf('t = %.4f sec.\n', curr_t);
            fprintf('\tphi_b = %f\n', obj.phi_b);
            fprintf('\tdV_ib_k = %f\n', obj.dV_ib_k);
            fprintf('\tgBias = %f\n', obj.gBias);
            fprintf('\taBias = %f\n', obj.aBias);
            fprintf('\n\n'); 
            
            obj.prev_t    = curr_t;
            obj.prev_acc  = modAccel;
            obj.prev_gyro = modGyro;
            obj.prev_vel  = newVel;
        end
        
        
        %% Main method where new GPS measurements (coords and velocity) are passed to the method
        %    coords = ( latitude (L(t)), longitude (lambda(t)) )', most probably should be in radians;
        %    time_stamp - sec., the time sample which the coordinates are obtained for;
        function obj = updateGPSData(obj, time_stamp, coords)
            if length(coords) == 2
                dLambda = coords(2) - obj.lambda_in + obj.w_ie * time_stamp;
                % 1. Find the RM^in_n(t) according to the Eq.(12)
                obj.RM_in_nt(1,1) = cos(dLambda);
                obj.RM_in_nt(1,2) = -sin(coords(1)) * sin(dLambda);
                obj.RM_in_nt(1,3) = -cos(coords(1)) * sin(dLambda);
                
                obj.RM_in_nt(2,1) = sin(obj.L_in) * sin(dLambda);
                obj.RM_in_nt(2,2) = cos(obj.L_in) * cos(coords(1)) + sin(obj.L_in) * sin(coords(1)) * cos(dLambda);
                obj.RM_in_nt(2,3) = cos(obj.L_in) * sin(coords(1)) - sin(obj.L_in) * cos(coords(1)) * cos(dLambda);
                
                obj.RM_in_nt(3,1) = -cos(obj.L_in) * sin(dLambda);
                obj.RM_in_nt(3,2) = sin(obj.L_in) * cos(coords(1)) - cos(obj.L_in) * sin(coords(1)) * cos(dLambda);
                obj.RM_in_nt(3,3) = sin(obj.L_in) * sin(coords(1)) + cos(obj.L_in) * cos(coords(1)) * cos(dLambda);
            else
                warning('IFA_Method.updateGPSData(): Wrong input coordinates vector dimension!');
            end
        end
        
        
        %% Setter for sampling interval value
        function setTS_IMU(obj, newTS)
            obj.TS_IMU = newTS;
        end      
        %% Getter for sampling interval value
        function output = getTS_IMU(obj)
            output = obj.TS_IMU;
        end
              
        %% Setter for sampling interval value
        function setTS_GPS(obj, newTS)
            obj.TS_GPS = newTS;
        end      
        %% Getter for sampling interval value
        function output = getTS_GPS(obj)
            output = obj.TS_GPS;
        end
        
        
        %% Setter for the carrier latitude value at the start time
        function set_L_in(obj, new_L_in)
            obj.L_in = new_L_in;
        end
        %% Getter for the carrier latitude value at the start time
        function output = get_L_in(obj)
            output = obj.L_in;
        end
        
        %% Setter for the carrier longitude value at the start time
        function set_lambda_in(obj, new_lambda_in)
            obj.lambda_in = new_lambda_in;
        end
        %% Getter for the carrier longitude value at the start time
        function output = get_lambda_in(obj)
            output = obj.lambda_in;
        end        
               
    end    
    
    
    methods (Access = private) 
        
        %% Method provides Skew-Symmetric matrix from the input vector
        function updateGravityVector(obj, newAcc)
            if length(newAcc) == 3
                obj.curr_grav = obj.alp_coef * obj.curr_grav + (1 - obj.alp_coef) * newAcc;
            else
                warning('IFA_Method.updateGravityVector(): Wrong input Acc vector dimension!');
            end
        end
        
        
        %% Method provides Skew-Symmetric matrix from the input vector
        function output = getSkewSym(obj, inpVect)
            output = zeros(3);
            if length(inpVect) == 3
                output = [      0     -inpVect(3)  inpVect(2);...
                           inpVect(3)      0      -inpVect(1);...
                           -inpVect(2) inpVect(1)      0         ];
            else
                warning('IFA_Method.getSkewSym(): Wrong input vector dimension!');
            end
        end
        
        %% Convertion from so(3) space to SO(3) space
        function output = getSO3_from_Euclid(obj, inp_euclid_vect)
            output = zeros(3);
            if (numel(inp_euclid_vect) == 3)
                %vect_phi = obj.so3_to_euclid(inp_euclid_vect);
                vect_phi = inp_euclid_vect;
                len_phi  = sqrt(vect_phi(1)^2 + vect_phi(2)^2 + vect_phi(3)^2);
                vect_p   = vect_phi / len_phi;
                
                output = cos(len_phi) * diag(ones(3,1)) + ...
                    (1 - cos(len_phi))* vect_p * (vect_p') + ...
                    sin(len_phi) * obj.getSkewSym(vect_p);

            else
                warning('IFA_Method.getSO3_from_Euclid(): Wrong input vector dimension!');
            end
        end
        
        %% Method provides the determination of the RM^b(t)_bi matrix based on the Eq.(11)
        function output = getCurrRM_bt(obj, currGyro)
            output = zeros(3);
            if length(currGyro) == 3                
                % The solution based on the solving the diff equation
                output = obj.RM_bt_ib * obj.getSO3_from_Euclid(-currGyro * obj.TS_IMU);
            else
                warning('IFA_Method.getCurrRM_bt(): Wrong input Gyro vector dimension!');
            end
        end
        
        
        %% Method provides Attitude determination vector in "ib" frame from the current Accel readings
        function output = getV_ib_k(obj, curr_t, currAccel, currGyro)
            output = zeros(3);
            if (length(currAccel) == 3) && (length(currGyro) == 3)
                dt = curr_t - obj.prev_t; 
                dw = currGyro - obj.prev_gyro;
                da = currAccel - obj.prev_acc;
                
                dTet1 = dt/2 * (dw/4 + currGyro - curr_t * dw/dt);
                dTet2 = dt*(dw/2 + currGyro) - curr_t * dw - dTet1;
                
                dV1 = dt/2 * (da/4 + currAccel - curr_t * da/dt);
                dV2 = dt*(da/2 + currAccel) - curr_t * da - dV1;
                
                dTet1_ss = obj.getSkewSym(dTet1);
                dTet2_ss = obj.getSkewSym(dTet2);
                dV1_ss   = obj.getSkewSym(dV1);
                
                output = dt/30 * (25*dV1 + 5*dV2 + 12*dTet1_ss*dV1 + 8*dTet1_ss*dV2 + ...
                    2*dV1_ss*dTet2 + 2*dTet2_ss*dV2);
            else
                warning('IFA_Method.getV_ib_k(): Wrong input Accel or Gyro vector dimension!');
            end
        end      
        
        
        %% Method provides Attitude determination vector in "in" frame from the current GPS readings
        function output = getV_in_k(obj, curr_t, gpsVelocity)
            output = zeros(3);
            if length(gpsVelocity) == 3                
                dt = curr_t - obj.prev_t;
                dV = gpsVelocity - obj.prev_vel;                
                w_n_ss = obj.getSkewSym(2*obj.w_n_ie + obj.w_n_en);
                
                term1 = eye(3) * gpsVelocity * dt + 0.5 * (dt^2) * (obj.w_n_in_ss * gpsVelocity + ...
                          eye(3) * dV / dt) + 1/3 * (dt^3) * obj.w_n_in_ss * dV / dt ;
                
                term2 = dt * obj.prev_vel;
                                
                term3 = (1/3*(dt^2)*eye(3) + 1/12*(dt^3)*obj.w_n_in_ss) * w_n_ss * obj.prev_vel + ...
                        (1/6*(dt^2)*eye(3) + 1/12*(dt^3)*obj.w_n_in_ss) * w_n_ss * gpsVelocity;
                
                term4 = 0.5 * (dt^2) * (eye(3) + dt / 3 * obj.w_n_in_ss) * obj.curr_grav;
                
                output = term1 + term2 + term3 + term4;
            else
                warning('IFA_Method.getV_in_k(): Wrong input GPS velocity or Navigation gravity vector dimension!');
            end
        end 
        
    end
end








