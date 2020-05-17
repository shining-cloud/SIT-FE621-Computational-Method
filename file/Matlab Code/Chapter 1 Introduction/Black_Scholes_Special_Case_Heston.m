% Illustration of Heston reducing to Black Scholes when sigma = 0

clc; clear;

% Strike Price, spot price, risk free rate, dividend yield, maturity
K = 100;
S = 100;
r = 0.03;
q = 0.02;
T = 0.5;

% Heston parameters 
% sigma and rho: not needed.  
% theta : Black Scholes variance
kappa = 5;
v0    = 0.05;
theta = v0;
lambda = 0;

% Integration grid
dphi = 0.01;
Uphi = 100;
Lphi = 0.00001;

%% The Black-Scholes Call and Put.  Uses v0 as the variance
d1 = (log(S/K) + (r-q+v0/2)*T)/sqrt(v0*T);
d2 = d1 - sqrt(v0*T);
BSCall = S*exp(-q*T)*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
BSPut = K*exp(-r*T)*normcdf(-d2) - S*exp(-q*T)*normcdf(-d1);

%% The Heston call and put
HCall = HestonPriceZeroSigma('C',kappa,theta,lambda,T,K,S,r,q,Uphi,dphi,Lphi);
HPut  = HestonPriceZeroSigma('P',kappa,theta,lambda,T,K,S,r,q,Uphi,dphi,Lphi);

%% Output the results
fprintf('Black Scholes as a special case of Heston\n')
fprintf('Using a volatility of %6.4f\n',theta);
fprintf('--------------------------------------------\n')
fprintf('Model               Call            Put \n');
fprintf('--------------------------------------------\n')
fprintf('Black Scholes    %10.8f %15.8f\n',BSCall,BSPut);
fprintf('Heston           %10.8f %15.8f\n',HCall,HPut);
fprintf('--------------------------------------------\n')

