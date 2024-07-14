% This script depicts raw ECG captured by the Arduido sensor
close all;
clear, clc;

%% ===================== PATHES SETTINGS ======================= %
% Date of the log-files
logFolderName = 'steady-pos';
logFileName   = 'MT_077002EC_000-000';

curDir = pwd;
pIncluder;

%% ====================== READING CONFIG ======================== %
% Read 'config.ini'
INI = INI('File','config.ini').read();


%% ======================= READING LOGS ========================= %
% Samples time - in seconds
fname = [logsFolder filesep logFileName '.' INI.general.logFilesExt];
%[timeIMU, accData, gyroData, anglesData] = LOGS_READER.readIMU(fname, INI);
[timeIMU, accData, gyroData] = LOGS_READER.readIMU(fname, INI);

fname = [logsFolder filesep logFileName '.' INI.general.logFilesExt];
[timeGPS, gpsLocData, gpsVelData] = LOGS_READER.readGPS(fname, INI);

% Calculate averaged IMU sample period of input data
TS = 0;
for t=2:length(timeIMU)
    TS = TS + timeIMU(t) - timeIMU(t-1);
end
dataParams.IMU_TS = TS / (length(timeIMU)-1);
dataParams.IMU_FS = 1 / dataParams.IMU_TS;

% Calculate averaged GPS sample period of input data
TS = 0;
for t=2:length(timeGPS)
    TS = TS + timeGPS(t) - timeGPS(t-1);
end
dataParams.GPS_TS = TS / (length(timeGPS)-1);
dataParams.GPS_FS = 1 / dataParams.GPS_TS;

%% 
if INI.debug.showDebugInfo
    fprintf('Averaged IMU sampling period: %5.3f sec.\n', dataParams.IMU_TS);
    fprintf('Averaged IMU sampling frequency: %5.3f Hz\n', dataParams.IMU_FS);
    fprintf('Averaged GPS sampling period: %5.3f sec.\n', dataParams.GPS_TS);
    fprintf('Averaged GPS sampling frequency: %5.3f Hz\n', dataParams.GPS_FS);
end
clear TS


%%
if INI.visualize.rawSensorPlots
    SENSORPLOTTER.plotOneSensorSignal(timeIMU, accData, 'Accelerometer', 'm/sec^2', ...
        {'acc_X', 'acc_Y', 'acc_Z'});
    SENSORPLOTTER.plotOneSensorSignal(timeIMU, gyroData, 'Gyroscope', 'rad/sec.', ...
        {'gyro_X', 'gyro_Y', 'gyro_Z'});
end

if INI.visualize.rawGPSDataPlots
    SENSORPLOTTER.plotOneSensorSignal(timeGPS, gpsLocData, 'GPS Location Data', 'degs.', ...
        {'Latitude', 'Longitude', 'Altitude'});
    SENSORPLOTTER.plotOneSensorSignal(timeGPS, gpsVelData, 'GPS Velocity Data', 'm/sec.', ...
        {'Vel_E', 'Vel_N', 'Vel_U'});
end


%% RM_b_n calculation method implementation
IFA_obj = IFA_Method('L_in', gpsLocData(1,1), 'lambda_in', gpsLocData(1,2),...
                     'TS_IMU', dataParams.IMU_TS, 'TS_GPS', dataParams.GPS_TS);
IFA_obj.Initialize(0, accData(1,:)', gyroData(1,:)', gpsVelData(1,:)');

time_GPS_index = 1;
gpsTsCheck = true;
for t=2:length(timeIMU)
    % If it is time for the GPS data, update the method with it
    if gpsTsCheck && (timeGPS(time_GPS_index) <= timeIMU(t))
        IFA_obj.updateGPSData(timeGPS(time_GPS_index), gpsLocData(time_GPS_index, :));
        time_GPS_index = time_GPS_index + 1;
        if time_GPS_index > length(time_GPS_index)
            gpsTsCheck = false;
        end
    end
    
    % Update the method with the new IMU data
    if time_GPS_index-1 == 0
        IFA_obj.updateIMUData(timeIMU(t), accData(t,:)', gyroData(t,:)', gpsVelData(1,:)');
    else
        IFA_obj.updateIMUData(timeIMU(t), accData(t,:)', gyroData(t,:)', gpsVelData(time_GPS_index-1,:)');
    end
    
    % IFA_obj.RM_b_n    % RM from B to N frames which calc-d for each timestamp
end





