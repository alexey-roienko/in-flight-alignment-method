%% State Extrapolation Equation
%% Output of the KF
 Theta_Fout_Kalman = zeros(length(time), 1);
 K_GainOutPut = zeros(length(time), 1);

%% Initialization values
 Xinit = 0;
 Pinit = 10;
 varQ  = 0.0003;
 varR  = 0.1;

%% Constant Matrices
 F = 1;
 B = TS;
 H = 1;
 stdG=std(gyroData) ; 
 stdA=std(accData) ; 
 Q= (stdG)*(stdG)'*TS'; 
 R= (stdA)*(stdA)'*TS'; 

%% Initialization Stage
 Theta_prv = Xinit;
 P = Pinit;

%% equation 1 
 for  t=1:1:length(time)
    U=gyroData(t,2);
   
    Theta_Pred=  F * Theta_prv+ B * U;
   %% Equation 2 
     
     P_Pred=F *P *F' + Q;    
   %% Equation 4a
     Z_k = accAngles(t,2);
     Y_bar= Z_k- H * Theta_Pred;
   %% Equation 3
    S = H * P_Pred * H' + R;
    K_Gain= H*P_Pred* H'/S;
 
   %% Equation 4
    Theta_Fout = Theta_Pred+ K_Gain* Y_bar;
   
   %% Equation 5
    I = eye(1);    % eye(size(X,1))
    P = (I-K_Gain)*P_Pred;
   
   %% Saving the current filter results
    Theta_Fout_Kalman(t, 1) = Theta_Fout ;
    K_GainOutPut(t, 1) =  K_Gain;  
    
 end
     
figure();
subplot(1,2,1);
plot(time, rad2deg(accAngles(:,2)), time, rad2deg( gyroAngles(:,2)),...
     time, rad2deg(Theta_Fout_Kalman), time, anglesData(:,2));
legend('Pitch Acc Angle', 'Pitch Gyro Angle', ' Pitch KF1', ' Pitch_GT angle');
grid on;

subplot(1,2,2);
plot(time,   K_GainOutPut);
legend('Kalman Gain')
grid on;
ylim([-0.05 1.05]);