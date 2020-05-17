% Illustration of pricing using uniform and non-uniform grids
clc; clear;

% Strike price, risk free rate, dividend yield, and maturity
K = 100;
r = 0.02;
q = 0.05;
Mat = 0.15;

% Heston parameters
kappa =  1.5;
theta =  0.04;
sigma =  0.3;
rho   = -0.9;
v0    =  0.05;
lambda = 0;
params = [kappa theta sigma v0 rho lambda];

% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0;  Smax = 2*K;
Vmin = 0;  Vmax = 0.5;
Tmin = 0;  Tmax = Mat;

% Number of grid points for the stock, volatility, and maturity
nS = 79;        % Stock price
nV = 39;        % Volatility
nT = 3000;      % Maturity

% The maturity time increment and grid
dt = (Tmax-Tmin)/nT;
T = [0:nT].*dt;

%% Pricing Using a Uniform Grid
% Increment for Stock Price and volatility
ds = (Smax-Smin)/nS;
dv = (Vmax-Vmin)/nV;

% The stock price and volatility grids
S = [0:nS].*ds;
V = [0:nV].*dv;

% Solve the PDE
U = HestonExplicitPDE(params,K,r,q,S,V,T);

% Obtain the price by 2-D interpolation
S0 = 101.52;
V0 = 0.05412;
UniformPrice = interp2(V,S,U,V0,S0);


%% Pricing Using a Non-Uniform Grid
% The stock price grid
c = K/5;
dz = 1/nS*(asinh((Smax-K)/c) - asinh(-K/c));
for i=1:nS+1;
	z(i) = asinh(-K/c) + (i-1)*dz;
	S(i) = K + c*sinh(z(i));
end

% The volatility grid
d = Vmax/500;
dn = asinh(Vmax/d)/nV;
for j=1:nV+1
	n(j) = (j-1)*dn;
	V(j) = d*sinh(n(j));
end

% Solve the PDE
U = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T);

% Obtain the price by 2-D interpolation
NonUniformPrice = interp2(V,S,U,V0,S0);


%% Closed form Price and errors
trap = 1;
PutCall = 'C';
[x w] = GenerateGaussLaguerre(32);
ClosedPrice = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);
UError = UniformPrice - ClosedPrice;
NError = NonUniformPrice  - ClosedPrice;


%% Output the results
clc;
fprintf('Stock price grid size  %5.0f\n', nS+1)
fprintf('Volatility grid size   %5.0f\n', nV+1)
fprintf('Number of time steps   %5.0f\n', nT)
fprintf('------------------------------------------\n')
fprintf('Method                Price   DollarError \n')
fprintf('------------------------------------------\n')
fprintf('Closed Form       %10.4f         \n', ClosedPrice)
fprintf('Uniform Grid      %10.4f    %5.2f\n', UniformPrice,UError)
fprintf('Non-Uniform Grid  %10.4f    %5.2f\n', NonUniformPrice,NError)
fprintf('------------------------------------------\n')
