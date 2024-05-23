classdef SENSORPLOTTER < handle
    %PLOTTER Summary of this class goes here
    %   Detailed explanation goes here      
    methods (Access = public, Static)
               
        %% Method plots three components of the single sensor signal defined as input argument
        function plotOneSensorSignal(time, signals, figName, yLabel, legendTitles, varargin)
            figure('units', 'normalized', 'outerposition', [0.05 0.05 .9 .9], 'IntegerHandle', 'off', ...
                'Name', figName);
            fSize = 14;            
            limits = zeros(4,2);
            
            limits(1,:) = [0 time(end,1)];
            
            % If Y limits are set by user...
            if ~isempty(varargin)
                limits(2,:) = varargin{1};
                limits(3,:) = varargin{1};
                limits(4,:) = varargin{1};
            end
                
            for s=1:3
                % Calculate Y limits automatically...
                if isempty(varargin)
                    minVal = min(signals(:,s));
                    maxVal = max(signals(:,s));
                    if maxVal > abs(minVal)
                        delta = 0.2*maxVal;
                        limits(s+1,:) = [-maxVal-delta maxVal+delta];
                    else
                        delta = 0.2*abs(minVal);
                        if delta < eps
                            limits(s+1,:) = [-0.1 0.1];
                        else
                            limits(s+1,:) = [minVal-delta abs(minVal)+delta];
                        end
                    end
                end
                
                subplot(3,1,s)
                plot(time(:, 1), signals(:, s));
                legend(legendTitles{s});
                xlabel('Time, sec.');
                ylabel(yLabel);
                xlim(limits(1,:));
                ylim(limits(s+1,:));
                set(gca,'fontsize',fSize);
                grid on;
                grid minor;
            end
        end
        
        
        %% Method plots three components of the two sensor signals defined as input arguments
        function plotTwoSensorSignals(time, signal1, signal2, figName, yLabel, legendTitles, varargin)
            figure('units', 'normalized', 'outerposition', [0.05 0.05 .9 .9], 'IntegerHandle', 'off', ...
                'Name', figName);
            fSize = 14;            
            limits = zeros(4,2);
            
            limits(1,:) = [0 time(end,1)];
            
            % If Y limits are set by user...
            if ~isempty(varargin)
                limits(2,:) = varargin{1};
                limits(3,:) = varargin{1};
                limits(4,:) = varargin{1};
            end
                
            for s=1:3
                % Calculate Y limits automatically...
                if isempty(varargin)
                    minVal1 = min(signal1(:,s));
                    minVal2 = min(signal2(:,s));
                    minVal  = min(minVal1, minVal2);
                    maxVal1 = max(signal1(:,s));
                    maxVal2 = max(signal2(:,s));
                    maxVal  = max(maxVal1, maxVal2);
                    if maxVal > abs(minVal)
                        delta = 0.2*maxVal;
                        limits(s+1,:) = [-maxVal-delta maxVal+delta];
                    else
                        delta = 0.2*abs(minVal);
                        limits(s+1,:) = [minVal-delta abs(minVal)+delta];
                    end
                end
                
                subplot(3,1,s)
                plot(time(:, 1), signal1(:, s), 'LineWidth', 1);
                hold on;
                plot(time(:, 1), signal2(:, s));
                legend(legendTitles{2*s-1}, legendTitles{2*s});
                xlabel('Time, sec.');
                ylabel(yLabel);
                xlim(limits(1,:));
                ylim(limits(s+1,:));
                set(gca,'fontsize',fSize);
                grid on;
                grid minor;
            end
        end
        
        
        %% Method plots three components of the three sensor signals defined as input arguments
        function plotThreeSensorSignals(time, signal1, signal2, signal3, figName, yLabel, legendTitles, varargin)
            figure('units', 'normalized', 'outerposition', [0.05 0.05 .9 .9], 'IntegerHandle', 'off', ...
                'Name', figName);
            fSize = 14;            
            limits = zeros(4,2);
            
            limits(1,:) = [0 time(end,1)];
            
            % If Y limits are set by user...
            if ~isempty(varargin)
                limits(2,:) = varargin{1};
                limits(3,:) = varargin{1};
                limits(4,:) = varargin{1};
            end
                
            for s=1:3
                % Calculate Y limits automatically...
                if isempty(varargin)
                    minVal1 = min(signal1(:,s));
                    minVal2 = min(signal2(:,s));
                    minVal3 = min(signal3(:,s));
                    minVal  = min([minVal1, minVal2, minVal3]);
                    
                    maxVal1 = max(signal1(:,s));
                    maxVal2 = max(signal2(:,s));
                    maxVal3 = max(signal2(:,s));
                    maxVal  = max([maxVal1, maxVal2, maxVal3]);
                    if maxVal > abs(minVal)
                        delta = 0.2*maxVal;
                        limits(s+1,:) = [-maxVal-delta maxVal+delta];
                    else
                        delta = 0.2*abs(minVal);
                        limits(s+1,:) = [minVal-delta abs(minVal)+delta];
                    end
                end
                
                subplot(3,1,s)
                plot(time(:, 1), signal1(:, s), time(:, 1), signal2(:, s), time(:, 1), signal3(:, s));
                legend(legendTitles{3*s-2}, legendTitles{3*s-1}, legendTitles{3*s});
                xlabel('Time, sec.');
                ylabel(yLabel);
                xlim(limits(1,:));
                ylim(limits(s+1,:));
                set(gca,'fontsize',fSize);
                grid on;
                grid minor;
            end
        end
        
        
        %% Method plots four components of the three sensor signals defined as input arguments
        function plotFourSensorSignals(time, signal1, signal2, signal3, signal4, figName, yLabel, legendTitles, varargin)
            figure('units', 'normalized', 'outerposition', [0.05 0.05 .9 .9], 'IntegerHandle', 'off', ...
                'Name', figName);
            fSize = 14;            
            limits = zeros(4,2);
            
            limits(1,:) = [0 time(end,1)];
            
            % If Y limits are set by user...
            if ~isempty(varargin)
                limits(2,:) = varargin{1};
                limits(3,:) = varargin{1};
                limits(4,:) = varargin{1};
            end
                
            for s=1:3
                % Calculate Y limits automatically...
                if isempty(varargin)
                    minVal1 = min(signal1(:,s));
                    minVal2 = min(signal2(:,s));
                    minVal3 = min(signal3(:,s));
                    minVal4 = min(signal4(:,s));
                    minVal  = min([minVal1, minVal2, minVal3, minVal4]);
                    
                    maxVal1 = max(signal1(:,s));
                    maxVal2 = max(signal2(:,s));
                    maxVal3 = max(signal3(:,s));
                    maxVal4 = max(signal4(:,s));
                    maxVal  = max([maxVal1, maxVal2, maxVal3, maxVal4]);
                    if maxVal > abs(minVal)
                        delta = 0.2*maxVal;
                        limits(s+1,:) = [-maxVal-delta maxVal+delta];
                    else
                        delta = 0.2*abs(minVal);
                        limits(s+1,:) = [minVal-delta abs(minVal)+delta];
                    end
                end
                
                subplot(3,1,s)
                plot(time(:, 1), signal1(:, s), time(:, 1), signal2(:, s), ...
                     time(:, 1), signal3(:, s), time(:, 1), signal4(:, s));
                legend(legendTitles{4*s-3}, legendTitles{4*s-2}, legendTitles{4*s-1}, legendTitles{4*s});
                xlabel('Time, sec.');
                ylabel(yLabel);
                xlim(limits(1,:));
                ylim(limits(s+1,:));
                set(gca,'fontsize',fSize);
                grid on;
                grid minor;
            end
        end
        
        
        %% Method plots 2D plot in squared window
        function plotSignalIn2D(Xcomp, Ycomp, figName, yLabel, legendTitles, varargin)
            figure('units', 'normalized', 'outerposition', [0.1 0.1 .45 .8], 'IntegerHandle', 'off', ...
                'Name', figName);
            fSize = 14;            
            limits = zeros(2,2);
            
            limits(1,:) = [0 Xcomp(end,1)];
            
            % If Y limits are set by user...
            if ~isempty(varargin)
                limits(2,:) = varargin{1};
            end
                
            % Calculate Y limits automatically...
            if isempty(varargin)
                minVal = min(Ycomp);
                maxVal = max(Ycomp);
                if maxVal > abs(minVal)
                    delta = 0.2*maxVal;
                    limits(2,:) = [-maxVal-delta maxVal+delta];
                else
                    delta = 0.2*abs(minVal);
                    if delta < eps
                        limits(2,:) = [-0.1 0.1];
                    else
                        limits(2,:) = [minVal-delta abs(minVal)+delta];
                    end
                end
            end

            plot(Xcomp, Ycomp);
            legend(legendTitles);
            xlabel('Time, sec.');
            ylabel(yLabel);
            xlim(limits(1,:));
            ylim(limits(2,:));
            set(gca,'fontsize',fSize);
            grid on;
            grid minor;
        end
        
    end
end
