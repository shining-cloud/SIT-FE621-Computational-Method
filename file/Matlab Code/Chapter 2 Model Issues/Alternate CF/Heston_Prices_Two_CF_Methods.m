% Demonstrate that prices are identical when
%   f2(phi) and f1(phi) are used, or
%   f2(phi) and f1(phi) = f2(phi-i) / (S*exp((r-q)T))

clc; clear;

T = 0.5;          % Time to maturity.
S = 100;          % Spot price.
r = 0.03;         % Risk free rate.
q = 0.02;         % Dividend yield
kappa  =  0.2;    % Heston parameter: mean reversion speed.
theta  =  0.05;   % Heston parameter: mean reversion level.
sigma  =  0.3;    % Heston parameter: volatility of vol
rho    = -0.8;    % Heston parameter: correlation
lambda =  0;      % Heston parameter: risk
v0     =  0.22;   % Heston parameter: initial variance.
PutCall = 'C';    % 'C'all or 'P'ut
trap = 1;         % 1 = "Little Trap" formulation
                  % 0 = Original Heston formulation

%% Generate prices using (f1,f2), and (f2,f2/S/exp[r-qT])
[x w] = GenerateGaussLaguerre(32);
K = [90:110];
for k=1:length(K)
	% f2 and f1
	CF = 1;
	Price1(k) = HestonPriceGaussLaguerre(PutCall,S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,CF);
	% f2 and f2/(S*exp((r-q)T))
	CF = 2;
	Price2(k) = HestonPriceGaussLaguerre(PutCall,S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,CF);
end

%% Output the results
fprintf('-------------------------------- \n')
fprintf('   Two c.f.         One c.f.     \n')
fprintf('   (f1,f2)    (f2,f2/S/exp[rqT]) \n')
fprintf('-------------------------------- \n')
for k=1:length(K)
    fprintf('%12.6f %15.6f \n',Price1(k),Price2(k));
end
fprintf('-------------------------------- \n')

