% Quadratic-Exponential and Moment-Matching schemes
% for the Heston model.

clc; clear;

% Option features
r = 0.03;           % Risk free rate
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

% Simulation features
NS = 5000;          % Number of stock price paths
NT = 100;           % Number of time steps per path

% Parameters for the QE algorithm
gamma1 = 0.5;
gamma2 = 0.5;
phic = 1.5;
MC = 1;
icdf = 2;

%% Quadratic Exponential and Moment Matching Simulation Prices
[S V QPrice] = QEPrice(params,PutCall,S0,K,Mat,r,q,NT,NS,gamma1,gamma2,MC,phic,icdf);

[S V MPrice] = MMPrice(params,PutCall,S0,K,Mat,r,q,NT,NS);

%% Exact price
Uphi = 50;
Lphi = 1e-5;
dphi = 0.1;
Trap = 1;
HPrice = HestonPriceConsol(PutCall,kappa,theta,lambda,rho,sigma,Mat,K,S0,r,q,v0,Trap,Lphi,Uphi,dphi);

%% Simulation Errors
Qerror = QPrice - HPrice;
Merror = MPrice - HPrice;

%% Display the results
fprintf('Simulation using %5.0f price paths and %4.0f time steps\n',NS,NT)
fprintf('---------------------------------------\n');
fprintf('Method               Price      Error \n');
fprintf('---------------------------------------\n');
fprintf('Exact            %10.4f         \n',HPrice);
fprintf('Quadratic-Exp    %10.4f  %10.4f \n',QPrice,Qerror);
fprintf('Moment-Matching  %10.4f  %10.4f \n',MPrice,Merror);
fprintf('---------------------------------------\n');


