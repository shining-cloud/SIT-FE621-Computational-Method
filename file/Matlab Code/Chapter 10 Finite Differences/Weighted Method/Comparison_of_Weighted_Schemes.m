% Comparison of various Weighted Scheme
clc; clear;

% Strike price, risk free rate, dividend yield, and maturity
K   = 100;
r   = 0.02;
q   = 0.05;
Mat = 0.15;

% Spot price and desired volatility at which to price the option
S0 = 101.52;
V0 = 0.05412;


%% Heston parameters.  Case 1 of Hout and Foulon (Table 1)
kappa =  1.5;
theta =  0.04;
sigma =  0.3;
rho   = -0.9;
lambda = 0;
params = [kappa theta sigma V0 rho lambda];

% Calculate the exact price
trap = 1;
PutCall = 'C';
[x w] = GenerateGaussLaguerre(32);
HPrice = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);

%% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0; Smax = 2.5*K;
Vmin = 0; Vmax = 0.5;
Tmin = 0; Tmax = Mat;

% Number of stock, vol, and time points
nS = 29;
nV = 29;
nT = 19;

% Choose the grid type
GridType = 'NonUniform';

if strcmp(GridType,'Uniform')
    % Increment for Stock Price, Volatility, and Maturity
    ds = (Smax-Smin)/nS;
    dv = (Vmax-Vmin)/nV;
    % Grid Vectors for the Stock Price, Volatility, and Maturity
    NS = nS+1;
    NV = nV+1;
    S = [0:NS-1].*ds;
    V = [0:NV-1].*dv;
else
    % The stock price grid
    c = K/5;
    dz = 1/nS*(asinh((Smax-K)/c) - asinh(-K/c));
    for i=1:nS+1;
        z(i) = asinh(-K/c) + (i-1)*dz;
        S(i) = K + c*sinh(z(i));
    end
    % The volatility grid
    d = Vmax/10;
    dn = asinh(Vmax/d)/nV;
    for j=1:nV+1
        n(j) = (j-1)*dn;
        V(j) = d*sinh(n(j));
    end
    NS = nS+1;
    NV = nV+1;
end

% Maturity grid -- common to both uniform and non-uniform grid types
NT = nT+1;
dt = (Tmax-Tmin)/nT;
T = [0:NT-1].*dt;

% Total number of S x V points
%N = NS*NV;

%% Weighted methods

% Obtain the "L" matrix
if strcmp(GridType,'Uniform')
    [derS derSS derV1 derV2 derVV derSV R] = BuildDerivatives(S,V,T);
else
    [derS derSS derV1 derV2 derVV derSV R] = BuildDerivativesNonUniform(S,V,T);
 end
L = (r-q).*derS + kappa.*theta.*derV1 - kappa.*derV2 + (1/2).*derSS + (1/2).*sigma^2*derVV + rho.*sigma.*derSV - r.*R;

% Identity matrix for the A and B matrices
I = eye(NS*NV);

% Explicit, Implicit, and Crank-Nicolson prices
thet = [0 1 0.5];
for k=1:3
    A = (I - thet(k).*dt*L);
    B = (I + (1-thet(k)).*dt.*L);
    if thet(k) == 0;
        invA = I;
    else
        invA = inv(A);
    end
    Price(k) = WeightedMethod(thet(k),S0,V0,K,S,V,T,A,invA,B);
    Error(k) = Price(k) - HPrice;
end


%% Output everything
disp([GridType ' grid type'])
fprintf('Stock price grid size of %5.0f \n', NS);
fprintf('Volatility grid size of  %5.0f \n', NV);
fprintf('Number of time steps     %5.0f \n', NT);
fprintf('-----------------------------------------------------\n');
fprintf('Method                     Price        Dollar Error \n');
fprintf('-----------------------------------------------------\n');
fprintf('Closed form          %12.4f \n',HPrice);
fprintf('Explicit method      %12.4f   %12.4f \n',Price(1),Error(1));
fprintf('Implicit method      %12.4f   %12.4f \n',Price(2),Error(2));
fprintf('Crank-Nicolson       %12.4f   %12.4f \n',Price(3),Error(3));
fprintf('-----------------------------------------------------\n');
 