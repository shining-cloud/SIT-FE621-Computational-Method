clc; clear;
clear;

%% Gauss Laguerre abscissas and weights
[xs ws] = GenerateGaussLaguerre(32);
[xt wt] = GenerateGaussLegendre(20);


%% Parameters and inputs as per Chiarella and Ziogas' paper, Table 1
S0 = 100;        % Spot price.
K  = 100;        % Strike
tau = .25;       % Time to maturity.
r = 0.01;        % Risk free rate.
q = 0.12;        % Dividend Yield
kappa  = 4;      % Heston parameter: mean reversion speed.
theta  = 0.09;   % Heston parameter: mean reversion level.
sigma  = 0.1;    % Heston parameter: volatility of vol
V0     = 0.04;   % Heston parameter: initial variance.
rho    = 0;      % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % "Little Trap" formulation
params = [kappa theta sigma V0 rho lambda];

%% Find b0 and b1 using Chiarella and Ziogas algorithm
EV = @(Vt,theta,kappa,T,t) (theta + (Vt - theta)*exp(-kappa*(T-t)));

% Starting values
Evt = EV(V0,theta,kappa,tau,0);
V00 = Evt + sigma/abs(kappa)*sqrt(kappa*theta/2);
V10 = Evt - sigma/abs(kappa)*sqrt(kappa*theta/2);

b00 = max(log(r*K/q),log(K));
b10 = 0;

% Tolerance and number of steps for Volterra-like algorithm
tol0 = 0.005;
tol1 = 0.005;
Ntau = 25;
Ntol = 1e-8;

% Double integration scheme for early exercise 
DoubleType = 'GLe';   % 'GLe' or 'Trapz'

% Number of grid points for the double trapezoidal rule
Nt = 15;

% Integration limits for maturity
a = 1e-10;
b = tau;

% Integration limits for stock price
c = 1e-10;
d = 100;

% Find the vectors b0(t) and b1(t)
[B0 B1] = findB(tau,params,K,r,q,V00,V10,b00,b10,xs,ws,xt,wt,Nt,Ntau,tol0,tol1,Ntol,a,b,c,d,DoubleType);
b0 = B0(end);
b1 = B1(end);

%% Find the prices
[AmerCZ EuroCZ]  = CZAmerCall(S0,tau,params,K,r,q,xs,ws,xt,wt,Nt,b0,b1,a,b,c,d,DoubleType);
 
% Find the European prices
PutCall = 'C';
trap = 1;
EuroClosed = HestonPriceGaussLaguerre(S0,K,tau,r,q,kappa,theta,sigma,lambda,V0,rho,PutCall,trap,xs,ws);

%% Obtain the American price using the explicit method
AmerPDE = 3.743150385654432;
EuroPDE = 3.530253898464514;
AmerCVPDE = EuroClosed + (AmerPDE - EuroPDE);


%% Output the results
fprintf('  b(0)    b(1)\n')
fprintf('--------------------\n')
fprintf('%7.4f %7.4f \n',b0,b1);
fprintf('------------------------------------------------------\n');
fprintf('                   European    American   AmericanCV\n');
fprintf('------------------------------------------------------\n');
fprintf('Closed           %10.4f               \n',EuroClosed);
fprintf('Explicit method  %10.4f %10.4f %10.4f \n',EuroPDE,AmerPDE,AmerCVPDE);
fprintf('Chiarella        %10.4f %10.4f        \n',EuroCZ ,AmerCZ);
fprintf('------------------------------------------------------\n');


%% Plot the double integrand
% funNum = 2;
% T = [0:.01:tau];
% X = [-75:75];
% [y Int] = DoubleTrapezoidal(params,S0,K,tau,r,q,b0,b1,X,T,funNum);
% 
% Int = Int(2:end,2:end);
% surf(Int)
% ylabel('Maturity')
% xlabel('Price')

%% Plot the boundary
dt = tau/Ntau;
X  = dt:dt:tau;
BX = exp(B0+B1.*V0);
plot(X,BX)
 

%% Obtain the American price using the explicit method

% % Minimum and maximum values for the Stock Price, Volatility, and Maturity
% Smin = 0;  Smax = 3*K;
% Vmin = 0;  Vmax = 0.5;
% Tmin = 0;  Tmax = tau;
% 
% % Number of grid points for the stock, volatility, and maturity
% nS = 64;        % Stock price
% nV = 34;        % Volatility
% nT = 5000;      % Maturity
% NS = nS+1;
% NV = nV+1;
% NT = nT+1;
% 
% % The maturity time increment and grid
% dt = (Tmax-Tmin)/nT;
% T = [0:NT-1].*dt;
% 
% % The stock price grid S(i)
% c = K/1;        % Value used by Int 'T Hout and Foulon
% dz = 1/(NS-1)*(asinh((Smax-K)/c) - asinh(-K/c));
% for i=1:NS;
%     z(i) = asinh(-K/c) + (i-1)*dz;
%     S(i) = K + c*sinh(z(i));
% end
% 
% % The volatility grid V(j)
% d = Vmax/10;   % Value used by Int 'T Hout and Foulon
% dn = asinh(Vmax/d)/(NV-1);
% for j=1:NV
%     n(j) = (j-1)*dn;
%     V(j) = d*sinh(n(j));
% end
% 
% % Obtain the PDE prices
% EuroAmer = 'A';
% U = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer);
% AmerPDE = interp2(V,S,U,V0,S0);
% 
% EuroAmer = 'E';
% U = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer);
% EuroPDE = interp2(V,S,U,V0,S0);
% 
% % PDE control variate price
% AmerCV = EuroClosed + (AmerPDE - EuroPDE);



