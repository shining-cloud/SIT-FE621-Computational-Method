% Obtaining the price of the Heston call or put

clc; clear;

% Option features
S = 100;         % Spot price
K = 100;         % Strike price
tau = .5;        % Maturity
r = 0.03;        % Risk free rate
q = 0.0;         % Dividend yield
kappa = 5;       % Heston parameter : reversion speed
sigma = 0.5;     % Heston parameter : volatility of variance
rho   = -0.8;    % Heston parameter : correlation
theta = 0.05;    % Heston parameter : reversion level
v0    = 0.05;    % Heston parameter : initial variance
lambda = 0;      % Heston parameter : risk preference
Trap = 0;        % 0 = Original Heston formulation
                 % 1 = Albrecher et al formulation

% Integration range				 
Lphi = 1e-10;     % Lower limit
dphi = 0.01;      % Increment
Uphi = 100;       % Upper limit

%% Obtain the Heston put and call using the Heston integrands
t = tic;
HPut   = HestonPrice('P',kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Trap,Lphi,Uphi,dphi);
HCall  = HestonPrice('C',kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Trap,Lphi,Uphi,dphi);
t1 = toc(t);

%% Obtain the Heston put and call using the consolidatedintegrands
t = tic;
HPutC  = HestonPriceConsol('P',kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Trap,Lphi,Uphi,dphi);
HCallC = HestonPriceConsol('C',kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Trap,Lphi,Uphi,dphi);
t2 = toc(t);

%% Output the result
clc;
fprintf('------------------------------------------------\n');
fprintf('Integrand       CallPrice   PutPrice    CalcTime\n');
fprintf('------------------------------------------------\n');
fprintf('Heston       %10.4f  %10.4f  %10.4f \n',HCall,HPut,t1);
fprintf('Consolidated %10.4f  %10.4f  %10.4f \n',HCallC,HPutC,t2);
fprintf('------------------------------------------------\n');

