clc; clear;
clear;

%% Gauss Laguerre abscissas and weights
[xs ws] = GenerateGaussLaguerre(32);
[xt wt] = GenerateGaussLegendre(32);

%% Parameters and inputs as per Chiarella and Ziogas' paper, Table 1
S = 100;         % Spot price.
K = 100;         % Strike
tau = .25;       % Time to maturity.
r = 0.01;        % Risk free rate.
q = 0.12;        % Dividend Yield
kappa  = 4;      % Heston parameter: mean reversion speed.
theta  = 0.09;   % Heston parameter: mean reversion level.
sigma  = 0.1;    % Heston parameter: volatility of vol
V0     = 0.04;   % Heston parameter: initial variance.
rho    = 0.0;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % "Little Trap" formulation
params = [kappa theta sigma V0 rho lambda];


%% Find the prices
% time integration limits
a = 1e-5;
b = tau;

% price integration limits
c = 1e-5;
d = 150;

% Number of points for the double trapezoidal rule
Nt = 150;

% b0 and b1 from the Chiarella and Ziogas algorithm
b0 = 4.723396421033418;
b1 = 0.173153592424732;

% Find the prices
[AmerGL EuroGL]  = CZAmerCall(S,tau,params,K,r,q,xs,ws,xt,wt,Nt,b0,b1,a,b,c,d,'GLe');
[AmerTZ EuroTZ]  = CZAmerCall(S,tau,params,K,r,q,xs,ws,xt,wt,Nt,b0,b1,a,b,c,d,'Trapz');

%% Output the results
fprintf('American call price using Chiarella and Ziogas method \n')
fprintf('----------------------------------------------------- \n')
fprintf('Double Integration Method       American    European  \n')
fprintf('----------------------------------------------------- \n')
fprintf('Gauss-Legendre      %20.5f %10.5f \n',AmerGL,EuroGL)
fprintf('Trapezoidal         %20.5f %10.5f \n',AmerTZ,EuroTZ)
fprintf('----------------------------------------------------- \n')


