% Multi domain integration in Zhu (2010)
clc; clear;

% Define the parameters and inputs
S = 100;         % Spot price.
K = 100;         % Strike
r = 0.05;        % Risk free rate.
q = 0.01;        % Dividend yield
PutCall = 'C';   % 'P'ut or 'C'all
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.05;    % Heston parameter: mean reversion level.
sigma = 0.3;     % Heston parameter: volatility of vol
lambda = 0;      % Heston parameter: risk.
v0 = .05;        % Heston parameter: initial variance.
rho = 0.45;      % Heston parameter: correlation
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
T = 1.5;         % Time to maturity.

% Gauss Legendre abscissas and weights
[xGLe wGLe] = GenerateGaussLegendre(32);

% The domain of integration and the tolerance
lo = 1e-10;
hi = 150;
N  = 10;
dA = (hi - lo)/N;
A = [lo:dA:hi];
tol = 1e-6;

%% Calculate the "true" option price using Newton-Cotes
N = 10000;
method = 3;
a = 1e-20;
b = 150;
tic
PriceSimpson = HestonPriceNewtonCoates(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);
tSimpson = toc;

%% Calculate the price using a multi-domain of integration
tic
[PriceMD Domain Npoints] = HestonPriceGaussLegendreMD(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLe,wGLe,A,tol);
tMD = toc;
errorMD = (PriceMD - PriceSimpson);

%% Calculate the price using Newton-Cotes formulas
a = Domain(1);
b = Domain(2);
tic
PriceNC = HestonPriceNewtonCoates(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,Npoints);
tNC = toc;
errorNC = (PriceNC - PriceSimpson);

%% Output the results
clc;
fprintf('Multi-Domain tolerance   %16.3e\n',tol);
fprintf('Lower integration limit  %16.3e\n',a);
fprintf('Upper integration limit  %12.3f\n',b);
fprintf('Number of integration points   %3.0f\n',Npoints);
fprintf('------------------------------------------------------------\n');
fprintf('Method                 Price      Error     ComputationTime\n');
fprintf('------------------------------------------------------------\n');
fprintf('%5.0f-point Simpson  %8.4f          %14.4f\n',N,PriceSimpson,tSimpson);
fprintf('Multi-Domain         %8.4f   %10.2e %10.4f\n',PriceMD,errorMD,tMD);
fprintf('Newton-Cotes         %8.4f   %10.2e %10.4f\n',PriceNC,errorNC,tNC);
fprintf('------------------------------------------------------------\n');

