% Heston (1993) Call price by Attari (2004) formula
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% By Fabrice Douglas Rouah

clc; clear;

% Define the parameters and inputs
S = 30;          % Spot price.
K = 20;          % Strike
T = 1/12;        % Time to maturity.
r = 0.01;        % Risk free rate.
q = 0.0;         % Dividend yield
kappa = 1.4;     % Heston parameter: mean reversion speed.
theta = 0.05;    % Heston parameter: mean reversion level.
sigma = 0.3;     % Heston parameter: volatility of vol
v0    = 0.05;    % Heston parameter: initial variance.
rho   = -0.8;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
PutCall = 'C';

% Abscissas and weights
[x w] = GenerateGaussLaguerre(32);

%% Price using Heston and Attari formulations
PriceHeston = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
PriceAttari = AttariPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);

% Display the results
fprintf('Method              Call Price\n');
fprintf('------------------------------------\n');
fprintf('Heston integrand     %10.8f \n',PriceHeston);
fprintf('Attari integrand     %10.8f \n',PriceAttari);
fprintf('------------------------------------\n');

%% Create and plot the integrands
phi = [1e-10:.001:20];
for k=1:length(phi);
	Attari(k) = AttariProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap);
	Heston(k) = HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
end
Zero = zeros(1,length(phi));

plot(phi,Attari,'k-',phi,Heston,'r-',phi,Zero,'k:')
legend('Attari Integrand', 'Heston Integrand')


