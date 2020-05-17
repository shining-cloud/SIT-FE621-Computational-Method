% Euler and Milstein discretization of the Heston model.

clc; clear;

% Option features
r = 0.03;         % Risk free rate
q = 0.02;           % Dividend yield
Mat = .25;          % Maturity in years
S0 = 100;           % Spot price
K  = 90;            % Strike price
PutCall = 'C';      % 'P'ut or 'C'all

% Heston parameters
kappa =  6.2;      % Variance reversion speed
theta =  0.06;     % Variance reversion level
sigma =  0.5;      % Volatility of Variance
rho   = -0.70;     % Correlation between Brownian motions
v0    =  0.03;     % Initial variance
lambda = 0;        % Risk parameter
params = [kappa theta sigma v0 rho lambda];

% Exact price using the "consolidated" form of the integrand
Uphi = 50;
Lphi = 1e-5;
dphi = 0.1;
Trap = 1;
HPrice = HestonPriceConsol(PutCall,kappa,theta,lambda,rho,sigma,Mat,K,S0,r,q,v0,Trap,Lphi,Uphi,dphi);

% Simulation features
N = 5000;          % Number of stock price paths
T = 100;           % Number of time steps per path
alpha = 0.5;       % Weight for explicit-implicit scheme
negvar = 'T';      % Use the truncation scheme for negative variances


%% Simulate the processes and obtain the option prices
schemeV = {'E' 'M' 'IM' 'WM'};
for k=1:length(schemeV)
	[S V F(k) Price(k)] = EulerMilsteinPrice(schemeV(k),negvar,params,PutCall,S0,K,Mat,r,q,T,N,alpha);
	error(k) = (HPrice-Price(k))/HPrice*100;
end

%% Display the results
fprintf('Simulation using %5.0f simulations and %4.0f time steps\n',N,T)
fprintf('----------------------------------------------------------------------\n')
fprintf('Method                       Price     PercentError  VarianceOverrides\n')
fprintf('----------------------------------------------------------------------\n')
fprintf('Exact with Gauss-Laguerre %10.4f\n',HPrice);
fprintf('Euler Scheme              %10.4f %10.4f %10.0f\n',Price(1),error(1),F(1));
fprintf('Milstein Scheme           %10.4f %10.4f %10.0f\n',Price(2),error(2),F(2));
fprintf('Implicit Milstein Scheme  %10.4f %10.4f %10.0f\n',Price(3),error(3),F(3));
fprintf('Weighted Implied Mlistein %10.4f %10.4f %10.0f\n',Price(4),error(4),F(4));
fprintf('----------------------------------------------------------------------\n')


