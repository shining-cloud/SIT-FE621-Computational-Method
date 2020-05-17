clc; clear;

% Market strikes
K = [120   125   130   135   140   145   150];
NK = length(K);

% Option settings
S0 = 137.14;
rf = 0.0010;
q  = 0.0068;
PutCall = 'C';

% RMSE parameters
kappa =  8.9931;
theta =  0.0571;
sigma =  2.0000;
v0    =  0.0405;
rho   = -0.7899;
lambda = 0;
trap = 1;

% Finely-spaced maturities
T = [45:5:360]./365;
NT = length(T);

% Gauss-Legendre abscissas and weights
[xGLe wGLe] = GenerateGaussLegendre(32);

% Multi-domain integration rule settings
lo = 1e-10;
hi = 500;
Ndomains  = 50 ;
dA = (hi - lo)/Ndomains;
A = [lo:dA:hi];
tol = 1e-8;

% Settings for the bisection algorithm
a = 0.01;
b = 3;
MaxIter = 1000;
Tol = 1e-5;

%% Generate the implied volatility and total variance
for t=1:NT
    for k=1:NK
        [Call(k,t) LoDomain HiDomain Npoints(k,t)] = HestonPriceGaussLegendreMD(PutCall,S0,K(k),T(t),rf,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLe,wGLe,A,tol);
        IV(k,t) = BisecBSIV(PutCall,S0,K(k),rf,q,T(t),a,b,Call(k,t),Tol,MaxIter);
        TotalVar(k,t) = IV(k,t)^2 * T(t);
    end
end

%% Plot the results
% Define the colormap "cool" and "autumn" also look nice
cc = hsv(1.5*NT);

% Plot the total variance
for t=1:NT
    plot(K,TotalVar(:,t),'color',cc(t,:))
    hold on
end
xlabel('Strike')
ylabel('Total Variance')
