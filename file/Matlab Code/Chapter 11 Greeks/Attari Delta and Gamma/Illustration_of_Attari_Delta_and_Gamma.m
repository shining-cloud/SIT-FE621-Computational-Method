% Delta and Gamma using Attari (2004) formula
% Compares finite differences with closed-form
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.

clc; clear;

% Define the parameters and inputs
S = 25;          % Spot price.
T = 1/2;         % Time to maturity.
r = 0.05;        % Risk free rate.
q = 0.0;         % Dividend yield
kappa = 1.4;     % Heston parameter: mean reversion speed.
theta = 0.05;    % Heston parameter: mean reversion level.
sigma = 0.5;     % Heston parameter: volatility of vol
v0    = .01;     % Heston parameter: initial variance.
rho   = -0.8;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
PutCall = 'C';

% Define the range of strikes
K = [20:30];

% Generate the abscissas and weights
[x w] = GenerateGaussLaguerre(32);


%% Finite difference Greeks
dS = 0.1;
for k=1:length(K);
    % Analytic Greeks
    Delta(k) = AttariGreeks(PutCall,kappa,theta,lambda,rho,sigma,T,K(k),S,r,q,v0,trap,'Delta',x,w);
    Gamma(k) = AttariGreeks(PutCall,kappa,theta,lambda,rho,sigma,T,K(k),S,r,q,v0,trap,'Gamma',x,w);
    % Greeks by finite differences
    Price1  = AttariPriceGaussLaguerre(PutCall,S+dS,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
    Price0  = AttariPriceGaussLaguerre(PutCall,S   ,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
    Price1_ = AttariPriceGaussLaguerre(PutCall,S-dS,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
    DeltaFD(k) = (Price1 - Price1_)/(2*dS);
    GammaFD(k) = (Price1 - 2*Price0 + Price1_)/dS^2;
end
Derror = sum(abs(Delta-DeltaFD));
Gerror = sum(abs(Gamma-GammaFD));


%% Output the results
fprintf('Delta and Gamma for the Attari (2004) model \n')
fprintf('-------------------------------------------------- \n');
fprintf('Strike   DeltaFD    Delta     GammaFD     Gamma    \n');
fprintf('-------------------------------------------------- \n');
for k=1:length(K);
    fprintf('%4.0f %10.4f %10.4f %10.4f %10.4f \n',K(k),DeltaFD(k),Delta(k),GammaFD(k),Gamma(k));
end
fprintf('-------------------------------------------------- \n');
fprintf('Delta sum absolute errors %8.2e \n',Derror)
fprintf('Gamma sum absolute errors %8.2e \n',Gerror)
fprintf('-------------------------------------------------- \n');

