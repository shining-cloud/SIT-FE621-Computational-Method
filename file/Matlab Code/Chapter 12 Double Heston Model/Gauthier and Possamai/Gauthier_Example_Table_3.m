% Double Heston example from P. Gauthier and D. Possamai
% "Efficient Simulation of the Double Heston Model"

clc; clear;

% Spot price, risk free rate, dividend yield
S = 61.90;
rf = 0.03;
q = 0;

% Parameters values
v01 = 0.6^2; 
v02 = 0.7^2;
sigma1 = 0.10;
sigma2 = 0.20;
kappa1 = 0.90;
kappa2 = 1.20;
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

% Define the strikes and maturities
percent = [1.0 1.0 0.7 0.7 1.3 1.3];
K = S.*percent;
T = [1 10 1 10 1 10];
PutCall = 'C';

%% Obtain call prices 
% Using the original Christoffersen, Heston, Jacobs characteristic function
% and using the "Little Trap" c.f.
a = 1e-20;
b = 100;
N = 500;
[x w] = GenerateGaussLaguerre(32);

for j=1:6
	OriginalTrapz(j) = DoubleHestonPriceTrapezoidal(PutCall,S,K(j),T(j),rf,q,param,1,a,b,N);
	GauthTrapz(j) = DoubleHestonPriceTrapezoidal(PutCall,S,K(j),T(j),rf,q,param,0,a,b,N);
	OriginalGLa(j) = DoubleHestonPriceGaussLaguerre(PutCall,S,K(j),T(j),rf,q,param,1,x,w);
	GauthGLa(j) = DoubleHestonPriceGaussLaguerre(PutCall,S,K(j),T(j),rf,q,param,0,x,w);
end

%% Display the call prices

fprintf('      Table 3 of Gauthier and Possamai (2010)\n')
fprintf('             Trapezoidal           Gauss-Laguerre\n')
fprintf('Strike   Original  LittleTrap   Original  LittleTrap\n')
fprintf('----------------------------------------------------\n')
for j=1:6
    fprintf('%5.2f  %10.4f %10.4f %10.4f %10.4f\n',K(j),OriginalTrapz(j),GauthTrapz(j),OriginalGLa(j),GauthGLa(j));
end
fprintf('----------------------------------------------------\n')

