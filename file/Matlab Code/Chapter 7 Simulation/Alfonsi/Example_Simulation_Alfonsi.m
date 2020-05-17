% Alfonsi discretization of the Heston model.

clc; clear;

% Option features
r = 0.03;           % Risk free rate
q = 0.02;           % Dividend yield
Mat = .25;          % Maturity in years
S0 = 100;           % Spot price
K  =  90;           % Strike price
PutCall = 'C';      % 'P'ut or 'C'all

% Heston parameters
kappa =  1.2;      % Variance reversion speed
theta =  0.06;     % Variance reversion level
sigma =  0.5;      % Volatility of Variance
rho   = -0.70;     % Correlation between Brownian motions
v0    =  0.03;     % Initial variance
lambda = 0;        % Risk parameter
trap = 1;          % "Little Trap" formulation
params = [kappa theta sigma v0 rho lambda];

% Exact price using the "consolidated" form of the integrand
[x w] = GenerateGaussLaguerre(32);
HPrice = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,params,trap,x,w);

% Simulation features
N = 5000;          % Number of stock price paths
T = 100;           % Number of time steps per path

%% Simulate the processes and obtain the option prices
[S V Price] = AlfonsiPrice(params,PutCall,S0,K,Mat,r,q,T,N);
DollarError = (HPrice-Price);

%% Display the results
fprintf('Simulation using %5.0f simulations and %4.0f time steps\n',N,T)
fprintf('-------------------------------------------------------\n')
fprintf('Method                       Price     DollarError  \n')
fprintf('-------------------------------------------------------\n')
fprintf('Exact with Gauss-Laguerre %10.4f\n',HPrice);
fprintf('Alfonsi Scheme            %10.4f %10.4f \n',Price,DollarError);
fprintf('-------------------------------------------------------\n')

