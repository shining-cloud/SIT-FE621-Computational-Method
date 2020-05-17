% Total volatility in the Heston model

clc; clear;

% Heston parameter settings
v0     =  0.4^2;     % Initial variance.
theta  =  0.2^2;     % Mean reversion level.
rho    = -0.6;       % Correlation
kappa  =  1.5;       % Mean reversion speed.
sigma  =  5.0;       % Volatility of variance.
lambda =  0;         % Risk.

S = 30;          % Spot price.
r = 0.0;         % Risk free rate.
q = 0.0;         % Dividend yield
PutCall = 'C';   % 'C'all or 'P'ut.
trap = 1;        % 1 = "Little Trap" form, 0 = Original Heston form

% Gauss Laguerre points
[x w] = GenerateGaussLaguerre(32);

% Generate the strikes and maturities.
K = [28:0.1:33];
T = [0.1:.05:1];


%% Generate the Heston prices, the implied volatilities, and the total volatilities.
a = 0.1;
b = 5;
Tol = 1e-5;
MaxIter = 1000;

for t=1:length(T);
    for k=1:length(K);
        Price = HestonPriceGaussLaguerre(PutCall,S,K(k),T(t),r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        IV(k,t) = BisecBSIV(PutCall,S,K(k),r,q,T(t),a,b,Price,Tol,MaxIter);
        TV(k,t) = IV(k,t) * sqrt(T(t));
    end
end

%% Surface plot of the Heston implied vol surface
Mat = round(T.*10)./10;
surf(IV)
box on
xlabel('Maturity')        
ylabel('Strike Price')
set(gca,'XTickLabel',Mat,'YTickLabel',K,'ZLim',[0.1 0.35]);


%% Plot the total volatility
%  plot(K,TV(:,1), K,TV(:,2), K,TV(:,3), K,TV(:,4), K,TV(:,5) , ...
%       K,TV(:,6), K,TV(:,7), K,TV(:,8), K,TV(:,9), K,TV(:,10), ...
%       K,TV(:,11),K,TV(:,12),K,TV(:,13),K,TV(:,14),K,TV(:,15), ...
%       K,TV(:,16),K,TV(:,17),K,TV(:,18),K,TV(:,9))
