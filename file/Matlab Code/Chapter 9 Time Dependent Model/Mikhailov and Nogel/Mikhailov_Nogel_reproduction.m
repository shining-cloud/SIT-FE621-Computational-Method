% Mikhailov and Nogel Example of Time Dependent Heston Model
% Reproduces Table 1 of their Wilmott magazine article
% Use Gauss-Laguerre integration

clc; clear;

% Option settings
S = 1;
rf = 0;
q = 0;
PutCall = 'C';

% Static parameters
theta  = 0.1;
sigma  = 0.2;
v0     = 0.1;
rho    = -0.3;
lambda = 0;

% Past parameter values for kappa=4 and kappa=2
param0 = [4 theta sigma v0 rho;...
	      2 theta sigma v0 rho];

% Current parameter values for kappa=1
param  = [1 theta sigma v0 rho];

% Assume the maturity of 5 yrs is divided into 3 equal maturity increments
tau0 = [5/3 5/3];

% Current maturity
tau  =  5/3;

% Weights and abscissas
[x w] = GenerateGaussLaguerre(32);

% Strikes
K = [0.5:.25:1.5];

%% Time dependent prices. Column 3 of Table 1 of Mihailov and Nogel
for k=1:5
	MN(k) = MNPriceGaussLaguerre(param,param0,tau,tau0,K(k),S,PutCall,rf,q,x,w);
end

% For comparison, the time-independent (constant parameter) prices
% Use average value of kappa and the maturity of 5 years
paramInd = [4 theta sigma v0 rho];
T = 5;
trap = 1;

for k=1:length(K)
	PriceInd(k) = MNPriceGaussLaguerre(paramInd,[],T,[],K(k),S,PutCall,rf,q,x,w);
end

% True values from the MN paper
TruePrice = [0.543017 0.385109 0.273303 0.195434 0.141210];

%% Output the results
fprintf('Reproduction of Table 1 of Mikhailov and Nogel (2003)\n')
fprintf('-------------------------------------------\n');
fprintf('Strike   TruePrice  TD-Price  Static Price \n');
fprintf('-------------------------------------------\n');
for k=1:5
    fprintf('%5.2f %10.4f %10.4f %10.4f \n',K(k),TruePrice(k),MN(k),PriceInd(k));
end
fprintf('-------------------------------------------\n');

errorMN  = TruePrice - MN;
errorInd = TruePrice - PriceInd;

%% Error plot
plot(K,errorMN,'rx-',K,errorInd,'kx-')
legend('MN error','Independent error')