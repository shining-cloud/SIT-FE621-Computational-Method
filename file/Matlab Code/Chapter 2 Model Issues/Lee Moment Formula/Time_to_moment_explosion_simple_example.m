clc; clear;

% Moment explosition settings
kappa = 0.1;
sigma = 0.3;
rho1 =  0.7;
rho2 =  0.0;
rho3 = -0.7;
lambda = 1;      % Moment explosition parameter
w = 2;           % Second moment

% Calculate the time of explosion of the second moment
T1 = MomentExplode(w,lambda,sigma,kappa,rho1);
T2 = MomentExplode(w,lambda,sigma,kappa,rho2);
T3 = MomentExplode(w,lambda,sigma,kappa,rho3);

%% Display the results
fprintf('Andersen and Piterbarg (2007) moment explosion times \n');
fprintf('------------------------------------------- \n');
fprintf('Correlation    Second Moment Explosion Time \n');
fprintf('------------------------------------------- \n');
fprintf('%8.2f   %20.2f \n',rho1,T1)
fprintf('%8.2f   %20.2f \n',rho2,T2)
fprintf('%8.2f   %20.2f \n',rho3,T3)
fprintf('------------------------------------------- \n');


