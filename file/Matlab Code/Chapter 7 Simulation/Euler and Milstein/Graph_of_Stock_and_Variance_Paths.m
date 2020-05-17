clc; clear;

% Option features
r = 0.03;           % Risk free rate
q = 0.02;           % Dividend yield
Mat = .25;          % Maturity in years
S0 = 100;           % Spot price
K  = 90;            % Strike price
PutCall = 'C';      % 'P'ut or 'C'all

% Heston parameters
kappa =  6;       % Variance reversion speed
theta =  0.06;    % Variance reversion level
sigma =  0.1;     % Volatility of Variance
v0    =  0.05;    % Initial variance
lambda = 0;       % Risk parameter

% Simulation features
N = 1;             % Number of stock price paths
T = 250;           % Number of time steps per path
alpha = 0.5;       % Weight for explicit-implicit scheme
negvar = 'T';      % Use the truncation scheme for negative variances
rho   =  0.9;      % Correlation between Brownian motions
params = [kappa theta sigma v0 rho lambda];


%% Simulate the processes and obtain the option prices
schemeV = 'M';
[S V F Price] = EulerMilsteinPrice(schemeV,negvar,params,PutCall,S0,K,Mat,r,q,T,N,alpha);

%% Plot the results
X = (1:T);
[a,h1,h2] = plotyy(X,S,X,V,'plot');
set(h1,'Color','k');
set(h2,'Color','r');
set(a(1),'YColor','k');
set(a(2),'YColor','r');
legend('Variance Level','Stock Price')
clear h1 h2 a
