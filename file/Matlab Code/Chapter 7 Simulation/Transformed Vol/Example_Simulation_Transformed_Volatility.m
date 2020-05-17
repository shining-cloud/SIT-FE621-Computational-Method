% Transformed Volatility discretization of the Heston model.

clc; clear;

% Option features
r = 0.03;           % Risk free rate
q = 0.02;           % Dividend yield
Mat = 3/12;         % Maturity in years
S0 = 100;           % Spot price
K  = 90;            % Strike price
PutCall = 'C';      % 'P'ut or 'C'all

% Heston parameters
kappa =  6.2;      % Variance reversion speed
theta =  0.06;     % Variance reversion level
sigma =  0.5;      % Volatility of Variance
rho   = -0.70;     % Correlation between Brownian motions
v0    =  0.03;     % Initial variance
sv0   = sqrt(v0);  % Initial volatility
lambda = 0;        % Risk parameter

%% Exact price using the "consolidated" form of the integrand
Uphi = 50;
Lphi = 1e-5;
dphi = 0.1;
Trap = 1;
params = [kappa theta sigma v0 rho lambda];
HPrice = HestonPriceConsol(PutCall,kappa,theta,lambda,rho,sigma,Mat,K,S0,r,q,v0,Trap,Lphi,Uphi,dphi);

%% Simulation features
NS = 5000;          % Number of stock price paths
NT = 100;           % Number of time steps per path
schemeV = {'Euler' 'TV'};

%% Run the simulations and display the results
sparams = [kappa theta sigma sv0 rho lambda];
for k=1:length(schemeV)
	[S v Price(k)] = TransVolPrice(schemeV(k),params,PutCall,S0,K,Mat,r,q,NT,NS);
	error(k) = HPrice-Price(k);
end


%% Display the results
clc;
fprintf('Simulation using %5.0f price paths and %4.0f time steps\n',NS,NT)
fprintf('---------------------------------------------------\n')
fprintf('Method                       Price       Error\n')
fprintf('---------------------------------------------------\n')
fprintf('Exact with Gauss-Laguerre %10.4f\n',HPrice);
fprintf('Euler Scheme              %10.4f %10.4f \n',Price(1),error(1));
fprintf('Transformed Vol Scheme    %10.4f %10.4f \n',Price(2),error(2));
fprintf('---------------------------------------------------\n')


