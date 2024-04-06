classdef FILTERS_TO < handle
    %% FILTERS_TO class includes different signal filters for raw IMU sensor data
    
    methods (Access = public, Static)
        
        %% Butterworth Low-pass filter wrapper
        function output = fButterworthFilter(input, INI)
            % Calculate Butterworth filter coefficients
            [b, a] = butter(INI.Accel.BtwFilterOrder,...
                    (2*INI.Accel.BtwFilterCutoffF) / (1/INI.General.dataSamplePeriod), 'low');
    
            % Apply the Butterworth filter
            output = zeros(size(input));
            output(:,1) = filter(b, a, input(:,1));
            output(:,2) = filter(b, a, input(:,2));
            output(:,3) = filter(b, a, input(:,3));
        end
        
        
        %% Low-pass Alpha-filter wrapper
        function output = fAlphaLPFilter(input, INI)
            % Apply the filter
            N = size(input,1);
            aCoef = INI.Accel.LPFilterAlpha;
            output = zeros(size(input));
            for a = 1:3
                output(1,a) = input(1,a);
                for s = 2:N
                    output(s,a) = aCoef * output(s-1,a) + (1-aCoef) * input(s,a);
                end
            end
        end
        
        
        %% Low-pass Alpha-filter wrapper
        function output = fComplimentaryFilter(gyroData, aAngles, filtCoef, TS)
            % Apply the filter
            N = size(aAngles,1);
            output = zeros(size(aAngles));
            output(1,:) = aAngles(1,:);
            for n = 2:N
                output(n,:) = filtCoef * output(n-1,:) + (1-filtCoef) * aAngles(n,:) + ...
                              filtCoef * gyroData(n,:) * TS;
            end
        end
    end
       
    
    
    methods (Access = private, Static)        
        %% Function implements averaging sliding window filtering
        % input   - data vector or matrix
        % WS      - sliding window size
        % modeStr - string described one of two modes:
        %           "round" - data filtering using rounding input samples
        %           "other str or empty" - ordinary filtering
        function output = aswFiltering(input, WS, modeStr)
            [X, Y] = size(input);
            output = zeros(X, Y);

            counter = WS;

            if ~(X==1 || Y==1)    
                fifo_vector = zeros(WS, Y);
            else
                N = length(input);
                fifo_vector = zeros(WS, 1);
            end

            if nargin==3
                if strcmpi(modeStr, 'round')
                    input = round(input);
                else                    
                    warning(['Undefined function parameter - ' modeStr]);
                    return
                end        
            end

            if ~(X==1 || Y==1)
                output(1:WS-1,:) = input(1:WS-1,:);
                fifo_vector(1:WS-1, :) = input(1:WS-1,:);
                for p=WS:X
                    if counter <= WS
                        fifo_vector(counter, :) = input(p,:);
                        counter = counter + 1;
                    else
                        counter = 1;
                        fifo_vector(counter, :) = input(p,:);
                    end    
                    output(p,:) = mean(fifo_vector,1);
                end
            else
                output(1:WS-1) = input(1:WS-1);
                fifo_vector(1:WS-1) = input(1:WS-1);
                for p=WS:N
                    if counter <= WS
                        fifo_vector(counter) = input(p);
                        counter = counter + 1;
                    else
                        counter = 1;
                        fifo_vector(counter) = input(p);
                    end    
                    output(p) = mean(fifo_vector);
                end
            end
        end
        
        
        %% Function implements averaging sliding window filtering
        % input   - data vector or matrix
        % WS      - sliding window size
        % modeStr - string described one of two modes:
        %           "round" - data filtering using rounding input samples
        %           "other str or empty" - ordinary filtering
        function output = fASWFiltering2(input, WS, modeStr)
            % Initialization
            [X, Y] = size(input);
            output = zeros(X, Y);
            fifoSize = 0;
            counter = 1;

            % Format definition
            if ~(X==1 || Y==1)    
                fifo_vector = zeros(WS, Y);
            else
                fifo_vector = zeros(WS, 1);
            end

            if nargin==3
                if strcmpi(modeStr, 'round')
                    input = round(input);
                else
                    warning(['Undefined function parameter - ' modeStr]);
                    return
                end
            end

            % Filtering
            if ~(X==1 || Y==1)
                for x = 1:X
                    if counter > WS
                        counter = 1;
                    end
                    % Angles transform from [-180; 180] plain into [0; 360] / [0; 360 plain]
                    if x > 1
                        for y = 1:Y
                            if input(x, y) - input(x-1, y) > 180
                                input(x, y) = input(x, y) - 360;
                            elseif input(x, y) - input(x-1, y) < -180
                                input(x, y) = input(x, y) + 360;
                            end
                        end
                    end

                    % Filling the buffer
                    fifo_vector(counter, :) = input(x, :);
                    for i = 1:length(fifo_vector(1, :))
                        for j = 1:length(fifo_vector(:, 1))
                            if fifo_vector(j, i) == 0
                                break;
                            end
                            fifoSize = j;
                        end
                    end
                    output(x, :) = mean(fifo_vector(1:fifoSize, :), 1);

                    % Inverse transform of angles
                    for i = 1:Y
                        if output(x, i) > 180
                            output(x, i) = -360 + output(x, i);
                        elseif output(x, i) < -180
                            output(x, i) = 360 + output(x, i);
                        end
                    end

                    % Increase counter
                    counter = counter + 1;
                end
            else
                for x = 1:X
                    if counter > WS
                        counter = 1;
                    end

                    % Angles transform from [-180; 180] plain into [0; 360] / [0; 360 plain]
                    if x > 1
                        if input(x) - input(x-1) > 180
                            input(x) = input(x) - 360;
                        elseif input(x) - input(x-1) < -180
                            input(x) = input(x) + 360;
                        end
                    end

                    % Filling the buffer
                    fifo_vector(counter) = input(x);
                    for i = 1:length(fifo_vector)
                        if fifo_vector(i) == 0
                            break;
                        end
                        fifoSize = i;
                    end
                    output(x) = mean(fifo_vector(1:fifoSize));

                    % Inverse transform of angles
                    if output(x) > 180
                        output(x) = -360 + output(x);
                    elseif output(x) < -180
                        output(x) = 360 + output(x);
                    end

                    % Increase counter
                    counter = counter + 1;
                end
            end
        end
    end
end

