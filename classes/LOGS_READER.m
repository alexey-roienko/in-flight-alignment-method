classdef LOGS_READER < handle
    
    methods (Access = public, Static)
                
        %% Public method for reading the settings from the INI-config file
        function INI = fReadSettings(fileSettings)
            INI = LOGS_READER.fReadINI(fileSettings);
        end
        
        
        %% Public method for reading different logs
        function [output, varargout] = read(file, varargin)
            [~, name, ext] = fileparts(file);
            ini = varargin{1};
            
            %% TactigonOne logs
            if strcmpi(ini.general.logsSource, 'TactigonOne')
                if contains(name, 'TO_') && strcmpi(ext, '.txt')
                    [output, varargout{1}, varargout{2}, varargout{3}, varargout{4}] = ...
                        LOGS_READER.fReadTOImuTXT(file);                
                else
                    error ('Read function cannot read the specified file. Please check your code.');
                end
            
            %% xSens logs
            elseif strcmpi(ini.general.logsSource, 'xSens')
                if contains(name, 'xSens_') && strcmpi(ext, '.txt')
                    [output, varargout{1}, varargout{2}, varargout{3}, varargout{4}] = ...
                        LOGS_READER.fReadXSensImuTXT(file);
                elseif contains(name, 'MT_') && strcmpi(ext, '.txt')
                    [output, varargout{1}] = LOGS_READER.fReadXSensGpsTXT(file);
                else
                    error ('Read function cannot read the specified file. Please check your code.');
                end
            end
        end
        
        %% Public method for reading different logs
        function [output, varargout] = readIMU(file, varargin)
            [~, name, ext] = fileparts(file);
            ini = varargin{1};
            
            %% xSens logs
            if strcmpi(ini.general.logsSource, 'xSens')
                if contains(name, 'MT_') && strcmpi(ext, '.txt')
                    [output, varargout{1}, varargout{2}] = ...
                        LOGS_READER.fReadXSensImuTXT(file);
                else
                    error ('Read function cannot read the specified file. Please check your code.');
                end
            end
        end
        
        
        %% Public method for reading different logs
        function [output, varargout] = readGPS(file, varargin)
            [~, name, ext] = fileparts(file);
            ini = varargin{1};
            
            %% xSens logs
            if strcmpi(ini.general.logsSource, 'xSens')
                if contains(name, 'MT_') && strcmpi(ext, '.txt')
                    [output, varargout{1}, varargout{2}] = LOGS_READER.fReadXSensGpsTXT(file);
                else
                    error ('Read function cannot read the specified file. Please check your code.');
                end
            end
        end
        
    end
    
    
    %% PRIVATE methods
    methods (Access = private, Static)
        
        %% Method reads settings from INI-file using open-source library
        function output = fReadINI(file)
            output = INI('File', file).read();        
        end   
        
        
        %% Method for reading xSens log with GPS and Velocity data
        function [time, varargout] = fReadXSensGpsTXT(filename, varargin)            
            %% Initialize variables
            delimiter = ',';
            endRow = inf;
            if nargin < 2
                startRow = 14;
            else
                startRow = varargin{1};
            end

            %% Open the text file.
            fileID = fopen(filename,'r');

            %% Reads header line
            for i=1:startRow
                headerLine = fgetl(fileID);
            end
            headerCells = textscan(headerLine, '%s', 'Delimiter', delimiter, ....
                                   'TextType', 'string', 'EmptyValue', NaN, ...
                                   'ReturnOnError', false, 'EndOfLine', '\r\n');
            clear headerLine
            
            % Look for UTC_time
            cellNmbs = zeros(4*3+1, 1);
            for i=1:length(headerCells{1})
                if strcmpi(headerCells{1}(i), 'UTC_Nano')
                    cellNmbs(1) = i;
                    break;
                end
            end
            
            %Look for Latitude etc
            stN = 1;
            if cellNmbs(1)~= 0,  stN = cellNmbs(1);  end
            for i=stN:length(headerCells{1})
                % GPS
                if strcmpi(headerCells{1}(i), 'Latitude'),  cellNmbs(2) = i;  end
                if strcmpi(headerCells{1}(i), 'Longitude'), cellNmbs(3) = i;  end
                if strcmpi(headerCells{1}(i), 'Altitude'),  cellNmbs(4) = i;  end
                % Velocity
                if strcmpi(headerCells{1}(i), 'Vel_E'),  cellNmbs(5) = i;  end
                if strcmpi(headerCells{1}(i), 'Vel_N'),  cellNmbs(6) = i;  end
                if strcmpi(headerCells{1}(i), 'Vel_U'),  cellNmbs(7) = i;  end
            end
            
            % Creating the 'formatSpec' using numbers found
            % If there are timesamples,...
            formatSpec = '';
            if cellNmbs(1) ~= 0
                for i = 1:cellNmbs(1)-1
                    formatSpec = [formatSpec '%s'];
                end
                % Add UTC_Nano values
                formatSpec = [formatSpec '%d'];
            end
            
            % Fill in gaps between parameters
            for i = cellNmbs(1)+1:cellNmbs(2)-1
                formatSpec = [formatSpec '%s'];
            end
            
            % Add GPS parameter values
            formatSpec = [formatSpec '%f%f%f'];

            % Look for Velocity values
            if (cellNmbs(4)+1 ~= cellNmbs(5))
                for i = cellNmbs(4)+1:cellNmbs(5)-1
                    formatSpec = [formatSpec '%s'];
                end
            end
            % Add Velocity values
            formatSpec = [formatSpec '%f%f%f'];

            % Add %s for remained values
            for i = cellNmbs(4)+1:length(headerCells{1})
                formatSpec = [formatSpec '%s'];
            end
            formatSpec = [formatSpec '%[^\n\r]'];

            
            %% Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow-startRow, 'Delimiter', delimiter, ....
                                 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', 0, ...
                                 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, ...
                    'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, ...
                    'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col}; dataArrayBlock{col}];
                end
            end

            %% Close the text file.
            fclose(fileID);

            %% Create output variable
            % If there is no UTC_Nano field, just output zeros for the 'time'
            if cellNmbs(1) == 0
                time = zeros(length(dataArray{1, cellNmbs(2)}), 1);
            else
                time = double(dataArray{1, cellNmbs(1)}) / 1e+9;
                
                time(1,1) = 0;
                for t=2:length(time)
                    % if the there is a new second, add 1 sec. to the current time sample
                    if (time(t-1) > time(t))                        
                        time(t:end) = 1 + time(t:end);
                    end                    
                end
            end
            
            % GPS data - Filtering the zero values
            temp = dataArray{1,cellNmbs(2)};
            nonzeroIndex = temp > 0;
            varargout{1}(:,1) = temp(nonzeroIndex);
            temp = dataArray{1,cellNmbs(3)};
            varargout{1}(:,2) = temp(nonzeroIndex);
            temp = dataArray{1,cellNmbs(4)};
            varargout{1}(:,3) = temp(nonzeroIndex);
            % Velocity data
            temp = dataArray{1,cellNmbs(5)};
            varargout{2}(:,1) = temp(nonzeroIndex);
            temp = dataArray{1,cellNmbs(6)};
            varargout{2}(:,2) = temp(nonzeroIndex);
            temp = dataArray{1,cellNmbs(7)};
            varargout{2}(:,3) = temp(nonzeroIndex);
            
            % Reduce the time samples
            time = time(nonzeroIndex);
        end
        
        
        %% Method for reading xSens log with Accel, Gyro, Magnet and
        %  Ground Truth Angles (from xSens board)        
        %       importfile(filename, startRow, endRow)
        function [time, varargout] = fReadXSensImuTXT(filename, varargin)            
            %% Initialize variables
            delimiter = ',';
            endRow = inf;
            if nargin<2
                startRow = 14;
            else
                startRow = varargin{1};
            end

            %% Open the text file.
            fileID = fopen(filename,'r');

            %% Reads header line
            for i=1:startRow
                headerLine = fgetl(fileID);
            end
            headerCells = textscan(headerLine, '%s', 'Delimiter', delimiter, ....
                                   'TextType', 'string', 'EmptyValue', NaN, ...
                                   'ReturnOnError', false, 'EndOfLine', '\r\n');
            clear headerLine
            
            % Look for UTC_time
            cellNmbs = zeros(4*3+1, 1);
            for i=1:length(headerCells{1})
                if strcmpi(headerCells{1}(i), 'UTC_Nano')
                    cellNmbs(1) = i;
                    break;
                end
            end
            
            %Look for Accel_X etc
            stN = 1;
            if cellNmbs(1)~= 0,  stN = cellNmbs(1);  end
            for i=stN:length(headerCells{1})
                % Accel
                if strcmpi(headerCells{1}(i), 'Acc_X'),  cellNmbs(2) = i;  end
                if strcmpi(headerCells{1}(i), 'Acc_Y'),  cellNmbs(3) = i;  end
                if strcmpi(headerCells{1}(i), 'Acc_Z'),  cellNmbs(4) = i;  end
                % Gyro
                if strcmpi(headerCells{1}(i), 'Gyr_X'),  cellNmbs(5) = i;  end
                if strcmpi(headerCells{1}(i), 'Gyr_Y'),  cellNmbs(6) = i;  end
                if strcmpi(headerCells{1}(i), 'Gyr_Z'),  cellNmbs(7) = i;  end
                % Magnet
                if strcmpi(headerCells{1}(i), 'Mag_X'),  cellNmbs(8) = i;  end
                if strcmpi(headerCells{1}(i), 'Mag_Y'),  cellNmbs(9) = i;  end
                if strcmpi(headerCells{1}(i), 'Mag_Z'),  cellNmbs(10) = i;  end
                % Euler angles
                if strcmpi(headerCells{1}(i), 'Roll'),   cellNmbs(11) = i;  end
                if strcmpi(headerCells{1}(i), 'Pitch'),  cellNmbs(12) = i;  end
                if strcmpi(headerCells{1}(i), 'Yaw'),    cellNmbs(13) = i;  end
            end
            
            % Creating the 'formatSpec' using numbers found
            % If there are timesamples,...
            formatSpec = '';
            if cellNmbs(1) ~= 0
                for i = 1:cellNmbs(1)-1
                    formatSpec = [formatSpec '%s'];
                end
                % Add UTC_Nano values
                formatSpec = [formatSpec '%d'];
            end
            
            % Step to the Accel values
            for i = cellNmbs(1)+1:cellNmbs(2)-1
                formatSpec = [formatSpec '%s'];
            end
            
            % Add Accel values
            formatSpec = [formatSpec '%f%f%f'];

            % Look for Gyro values
            if (cellNmbs(4)+1 ~= cellNmbs(5))
                for i = cellNmbs(4)+1:cellNmbs(5)-1
                    formatSpec = [formatSpec '%s'];
                end
            end
            % Add Gyro values
            formatSpec = [formatSpec '%f%f%f'];

            % Look for Magnet values
            if (cellNmbs(7)+1 ~= cellNmbs(8))
                for i = cellNmbs(7)+1:cellNmbs(8)-1
                    formatSpec = [formatSpec '%s'];
                end
            end
            % Add Magnet values
            formatSpec = [formatSpec '%f%f%f'];

            % Look for Euler angle values
            if (cellNmbs(10)+1 ~= cellNmbs(11))
                for i = cellNmbs(10)+1:cellNmbs(11)-1
                    formatSpec = [formatSpec '%s'];
                end
            end
            % Add Euler angle values
            formatSpec = [formatSpec '%f%f%f'];

            % Add %s for remained values
            for i = cellNmbs(13)+1:length(headerCells{1})
                formatSpec = [formatSpec '%s'];
            end
            formatSpec = [formatSpec '%[^\n\r]'];

            
            %% Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow-startRow, 'Delimiter', delimiter, ....
                                 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', 0, ...
                                 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, ...
                    'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, ...
                    'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col}; dataArrayBlock{col}];
                end
            end

            %% Close the text file.
            fclose(fileID);

            %% Create output variable
            % If there is no UTC_Nano field, just output zeros for the 'time'
            if cellNmbs(1) == 0
                time = zeros(length(dataArray{1, cellNmbs(2)}), 1);
            else
                time = double(dataArray{1, cellNmbs(1)}) / 1e+9;
                
                time(1,1) = 0;
                for t=2:length(time)
                    % if the there is a new second, add 1 sec. to the current time sample
                    if (time(t-1) > time(t))                        
                        time(t:end) = 1 + time(t:end);
                    end                    
                end
            end
            
            % Accel data
            varargout{1}(:,1) = dataArray{1,cellNmbs(2)};
            varargout{1}(:,2) = dataArray{1,cellNmbs(3)};
            varargout{1}(:,3) = dataArray{1,cellNmbs(4)};
            % Gyro data
            varargout{2}(:,1) = dataArray{1,cellNmbs(5)};
            varargout{2}(:,2) = dataArray{1,cellNmbs(6)};
            varargout{2}(:,3) = dataArray{1,cellNmbs(7)};
%             % Magnet data
%             varargout{3}(:,1) = dataArray{1,cellNmbs(8)};
%             varargout{3}(:,2) = dataArray{1,cellNmbs(9)};
%             varargout{3}(:,3) = dataArray{1,cellNmbs(10)};
%             % Angles data
%             varargout{4}(:,1) = dataArray{1,cellNmbs(11)};
%             varargout{4}(:,2) = dataArray{1,cellNmbs(12)};
%             varargout{4}(:,3) = dataArray{1,cellNmbs(13)};
        end

        
        
        %% Method for reading Tactigon One log with Accel, Gyro, Magnet and
        %  Ground Truth Angles from Tactigon One board
        function [time, varargout] = fReadTOImuTXT(fname)            
            % Initialize variables.
            if nargin<=2
                startRow = 1;
                endRow = inf;
            end

            % Format for each line of text:
            %   column6: double (%f)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%*6s%6f%*10s%6f%6f%6f%*10s%6f%6f%6f%*10s%5f%5f%5f%*12s%5f%5f%5f%[^\n\r]';

            % Open the text file.
            fileID = fopen(fname,'r');

            % Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '',...
                'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, ...
                'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1,...
                    'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN,...
                    'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                dataArray{1} = [dataArray{1}; dataArrayBlock{1}];
            end

            % Close the text file.
            fclose(fileID);

            % Post processing for unimportable data.
            % No unimportable data rules were applied during the import, so no post
            % processing code is included. To generate code which works for
            % unimportable data, select unimportable cells in a file and regenerate the
            % script.

            % Create output variable
            time(1,1) = 0;
            minCounter = 0;
            for t=2:length(dataArray{1,1}(:))
                % if the there is a new minute, add 60 sec. to the current
                % time sample
                if (dataArray{1,1}(t-1) > dataArray{1,1}(t))
                    minCounter = minCounter + 1;
                end
                time(t,1) = (60*minCounter + dataArray{1,1}(t)) - dataArray{1,1}(1);
            end
            % Accel data
            varargout{1}(:,1) = dataArray{1,2};
            varargout{1}(:,2) = dataArray{1,3};
            varargout{1}(:,3) = dataArray{1,4};       
            % Gyro data
            varargout{2}(:,1) = dataArray{1,5};
            varargout{2}(:,2) = dataArray{1,6};
            varargout{2}(:,3) = dataArray{1,7};       
            % Magnet data
            varargout{3}(:,1) = dataArray{1,8};
            varargout{3}(:,2) = dataArray{1,9};
            varargout{3}(:,3) = dataArray{1,10};            
            % Angles data
            varargout{4}(:,1) = dataArray{1,11};
            varargout{4}(:,2) = dataArray{1,12};
            varargout{4}(:,3) = dataArray{1,13};
        end
        
        
        
        %% Method for reading Tactigon One Accelerometer log
        function [time, varargout] = fReadTOAccelTXT(fname)            
            % Initialize variables.
            if nargin<=2
                startRow = 1;
                endRow = inf;
            end

            % Format for each line of text:
            %   column6: double (%f)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%*6s%6f%*13s%5f%*6s%5f%*6s%5f%[^\n\r]';

            % Open the text file.
            fileID = fopen(fname,'r');

            % Read columns of data according to the format.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '',...
                'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, ...
                'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1,...
                    'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN,...
                    'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                dataArray{1} = [dataArray{1}; dataArrayBlock{1}];
            end

            % Close the text file.
            fclose(fileID);

            % Post processing for unimportable data.
            % No unimportable data rules were applied during the import, so no post
            % processing code is included. To generate code which works for
            % unimportable data, select unimportable cells in a file and regenerate the
            % script.

            % Create output variable
            time(1,1) = 0;
            for t=2:length(dataArray{1,1}(:))
                % if the there is a new minute, add 60 sec. to the current
                % time sample
                if (dataArray{1,1}(t-1) > dataArray{1,1}(t))
                    time(t,1) = (60 + dataArray{1,1}(t)) - dataArray{1,1}(1);
                else
                    time(t,1) = dataArray{1,1}(t) - dataArray{1,1}(1);
                end
            end            
            varargout{1}(:,1) = dataArray{1,2};
            varargout{1}(:,2) = dataArray{1,3};
            varargout{1}(:,3) = dataArray{1,4};            
        end
        
        
        %% Method for reading an image file
        function im = fReadImage(path, imName)
            imPath = [path '\' imName];
            if (~isnan(imPath))
                im  = imread(imPath);
            end
        end
        
    end
end