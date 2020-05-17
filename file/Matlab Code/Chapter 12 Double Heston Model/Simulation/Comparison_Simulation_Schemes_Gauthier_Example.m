% Double Heston example from P. Gauthier and D. Possamai
% "Efficient Simulation of the Double Heston Model"

clc; clear;

% Spot price, risk free rate, dividend yield
S0 = 61.90;
K  = 61.90;
Mat = 1;
rf = 0.03;
q  = 0;
PutCall = 'C';

% Double Heston parameter values
v01 = 0.6^2; 
v02 = 0.7^2;
sigma1 = 0.1;
sigma2 = 0.2;
kappa1 = 0.9;
kappa2 = 1.2;
rho1 = -0.5;
rho2 = -0.5;
theta1 = 0.10;
theta2 = 0.15;

% Stack the parameters into a single vector
param(1) = kappa1;
param(2) = theta1;
param(3) = sigma1;
param(4) =    v01;
param(5) =   rho1;
param(6) = kappa2;
param(7) = theta2;
param(8) = sigma2;
param(9) =    v02;
param(10)=   rho2;

% Simulation settings
NS = 5000;
NT = 100;


%% Obtain true price and simulated prices
trap = 1;
[x w] = GenerateGaussLaguerre(32);
True = DoubleHestonPriceGaussLaguerre(PutCall,S0,K,Mat,rf,q,param,x,w,trap);

tic;
[S V1 V2 Euler] = DHEulerAlfonsiSim('Euler',param,S0,K,Mat,rf,q,NT,NS,PutCall);
t2=toc;
Eerror = True - Euler;
Ep = Eerror/True*100;
fprintf('Finished Euler Scheme in      %6.4f seconds\n', t2)

tic;
[S V1 V2 Alfonsi] = DHEulerAlfonsiSim('Alfonsi',param,S0,K,Mat,rf,q,NT,NS,PutCall);
t1=toc;
Aerror = True - Alfonsi;
Ap = Aerror/True*100;
fprintf('Finished Alfonsi Scheme in    %6.4f seconds\n', t1)


tic;
[S v1 v2 ZEuler] = DHTransVolSim('ZhuEuler',param,S0,K,Mat,rf,q,NT,NS,PutCall);
t3=toc;
Zerror = True - ZEuler;
Zp = Zerror/True*100;
fprintf('Finished Zhu-Euler scheme in  %6.4f seconds\n', t3)

tic;
[S v1 v2 TV] = DHTransVolSim('ZhuTV',param,S0,K,Mat,rf,q,NT,NS,PutCall);
t4=toc;
Terror = True - TV;
Tp = Terror/True*100;
fprintf('Finished Zhu-T.V. scheme in   %6.4f seconds\n', t4)

tic;
[S V1 V2 QE] = DHQuadExpSim(param,S0,K,Mat,rf,q,NT,NS,PutCall);
t5=toc;
Qerror = True - QE;
Qp = Qerror/True*100;
fprintf('Finished Quad Exp scheme in   %6.4f seconds\n', t5)

%% Display the results
fprintf('------------------------------------------------------------\n')
fprintf('Double Heston simulation schemes\n')
fprintf('Number of simuations: %5.0f\n',NS);
fprintf('Number of time steps: %5.0f\n',NT);
fprintf('------------------------------------------------------------\n')
fprintf('Method            Price   DollarError  PercentError  SimTime\n')
fprintf('------------------------------------------------------------\n')
fprintf('Closed Form       %6.4f                       \n', True)
fprintf('Euler             %6.4f  %8.3f  %10.3f %10.3f \n', Euler,Eerror,Ep,t2)
fprintf('Alfonsi           %6.4f  %8.3f  %10.3f %10.3f \n', Alfonsi,Aerror,Ap,t1)
fprintf('Zhu Euler         %6.4f  %8.3f  %10.3f %10.3f \n', ZEuler,Zerror,Zp,t3)
fprintf('Zhu Trans Vol     %6.4f  %8.3f  %10.3f %10.3f \n', TV,Terror,Tp,t4)
fprintf('Quadratric Exp    %6.4f  %8.3f  %10.3f %10.3f \n', QE,Qerror,Qp,t5)
fprintf('------------------------------------------------------------\n')


