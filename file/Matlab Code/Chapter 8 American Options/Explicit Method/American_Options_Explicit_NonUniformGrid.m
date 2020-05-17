% Solving the PDE for the Heston model for an American or European Call or Put.
% In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
% in the Heston Model with Correlation" Int J of Num Analysis and
% Modeling, 2010.
% Uses UNEVEN grid sizes

clc; clear;

% Settings from Clarke and Parrott
Spot = [8 9 10 11 12];
TruePrice = [2.00 1.107641 0.520030 0.213668 0.082036];
K = 10;
r = 0.1;
q = 0;
Mat = 0.25;
kappa = 5;
theta = 0.16;
sigma = 0.9;
rho = 0.1;
v0 = 0.0625;
V0 = v0;
lambda = 0;

% Option flavor
PutCall = 'P';
EuroAmer = 'A';

% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0;  Smax = 3*K;
Vmin = 0;  Vmax = 0.5;
Tmin = 0;  Tmax = Mat;

% Number of grid points for the stock, volatility, and maturity
nS = 64;        % Stock price
nV = 34;        % Volatility
nT = 5000;      % Maturity

NS = nS+1;
NV = nV+1;
NT = nT+1;

% The maturity time increment and grid
dt = (Tmax-Tmin)/nT;
T = [0:NT-1].*dt;


%% Loop through and find the PDE prices and the error

% The stock price grid S(i)
c = K/1;        % Value used by Int 'T Hout and Foulon
dz = 1/(NS-1)*(asinh((Smax-K)/c) - asinh(-K/c));
for i=1:NS;
    z(i) = asinh(-K/c) + (i-1)*dz;
    S(i) = K + c*sinh(z(i));
end

% The volatility grid V(j)
d = Vmax/10;   % Value used by Int 'T Hout and Foulon
dn = asinh(Vmax/d)/(NV-1);
for j=1:NV
    n(j) = (j-1)*dn;
    V(j) = d*sinh(n(j));
end

% Solve the PDE using the grid
params = [kappa theta sigma v0 rho lambda];
tic
U = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer);
toc

%% Calculate the prices from the PDE grid
disp('Grid Sizes')
fprintf('Stock price %3.0f steps, Volatility %3.0f steps, Time %5.0f steps \n',NS,NV,NT);
fprintf('---------------------------------------------------------------------\n')
fprintf('Spot Price  Clarke-Parrott    PDE      Dollar Error  Percentage Error\n')
fprintf('---------------------------------------------------------------------\n')
for k=1:length(Spot);
    S0 = Spot(k);
    AmerPrice(k) = interp2(V,S,U,V0,S0);
    DolError(k) = (TruePrice(k) - AmerPrice(k));
    RelError(k) = DolError(k)/TruePrice(k)*100;
    fprintf('%5.0f  %15.4f   %10.4f   %10.4f   %10.4f\n',Spot(k),TruePrice(k),AmerPrice(k),DolError(k),RelError(k));
end
fprintf('---------------------------------------------------------------------\n')



