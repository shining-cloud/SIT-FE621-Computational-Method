% Kahl-Jackel discretization of the Heston model.
% Compares the IJK, Pathwise simulation, and Balanced Implicit Schemes

%clc; clear;
clear;

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
N = 5000;          % Number of stock price paths
T = 100;           % Number of time steps per path
negvar = 'R';      % Negative variance scheme

%% IJK Price
scheme = 'IJK';
[S V Fi IPrice] = KahlJackelPrice(scheme,negvar,params,PutCall,S0,K,Mat,r,q,T,N);

% Pathwise Adapdated Linearization Quadratic Price
scheme = 'PW';
[S V Fp PPrice] = KahlJackelPrice(scheme,negvar,params,PutCall,S0,K,Mat,r,q,T,N);

% Balanced Implicit Price
scheme = 'B';
[S V Fb BPrice] = KahlJackelPrice(scheme,negvar,params,PutCall,S0,K,Mat,r,q,T,N);

% Exact price
Uphi = 50;
Lphi = 1e-5;
dphi = 0.1;
Trap = 1;
HPrice = HestonPriceConsol(PutCall,kappa,theta,lambda,rho,sigma,Mat,K,S0,r,q,v0,Trap,Lphi,Uphi,dphi);

%% Simulation Errors
Ierror = (IPrice-HPrice)/HPrice*100;
Perror = (PPrice-HPrice)/HPrice*100;
Berror = (BPrice-HPrice)/HPrice*100;


%% Display the results
clc;
fprintf('Simulation using %5.0f simulations and %4.0f time steps\n',N,T)
fprintf('-----------------------------------------------------\n')
fprintf('Scheme        Price      Error    # Negative Variances\n')
fprintf('-----------------------------------------------------\n')
fprintf('Exact      %10.4f\n',HPrice);
fprintf('IJK        %10.4f  %10.4f %10.0f \n', IPrice,Ierror,Fi)
fprintf('Pathwise   %10.4f  %10.4f %10.0f \n', PPrice,Perror,Fp)
fprintf('Balanced   %10.4f  %10.4f %10.0f \n', BPrice,Berror,Fb)
fprintf('-----------------------------------------------------\n')


