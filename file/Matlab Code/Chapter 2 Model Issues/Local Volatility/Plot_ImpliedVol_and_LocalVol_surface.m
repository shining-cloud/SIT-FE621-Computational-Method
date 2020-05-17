clc; clear;

% Heston parameter settings
v0    = 0.4^2;     % Initial variance.
theta = 0.2^2;     % Mean reversion level.
rho   = -0.6;      % Correlation
kappa = 1.5;       % Mean reversion speed.
sigma = 5.1;       % Volatility of variance.
lambda = 0;        % Risk.

% Option settings
S = 100;         % Spot price.
r = 0.0;         % Risk free rate.
q = 0.0;         % Dividend yield
PutCall = 'C';   % 'C'all or 'P'ut.
trap = 1;        % 1 = "Little trap" form, 0 = Original Heston form

% Gauss Laguerre points
[x w] = GenerateGaussLaguerre(32);

% Generate the strikes and maturities.
K = [95:0.5:105];
T = [0.25:.04:1.03];
NK = length(K);
NT = length(T);

%% Generate the Heston prices and the implied and local volatilities
a = 0.1;
b = 5;
Tol = 1e-5;
MaxIter = 1000;
dt = 0.001;
dK = .0125;

for t=1:NT;
    for k=1:NK;
        % Heston prices
        Price = HestonCallGaussLaguerre(S,K(k),T(t),r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        % Implied volatility
        IV(k,t) = BisecBSIV(PutCall,S,K(k),r,T(t),a,b,Price,Tol,MaxIter);
        % Local variance
        CT  = HestonCallGaussLaguerre(S,K(k),T(t)+dt,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        CT_ = HestonCallGaussLaguerre(S,K(k),T(t)-dt,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        dCdT =  (CT - CT_) / (2*dt);
        CK  = HestonCallGaussLaguerre(S,K(k)+dK,T(t),r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        CK0 = HestonCallGaussLaguerre(S,K(k)   ,T(t),r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        CK_ = HestonCallGaussLaguerre(S,K(k)-dK,T(t),r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        dC2dK2 = (CK - 2*CK0 + CK_) / (dK)^2;
        LocalVar = 2*dCdT / (K(k)^2*dC2dK2);
        % Local volatility
        LV(k,t) = sqrt(LocalVar);
    end
end

%% Surface plot of the Heston implied vol and local vol surface

% Transparent mesh for implied volatility
mesh(IV)
axis tight
view([100,150,100])

box on
xlabel('Maturity')        
ylabel('Strike Price')
set(gca,'XTick',      [0  5 10 15 20]);      % Maturity
set(gca,'XTickLabel', [.2  .4  .6  .8  1 ])        
set(gca,'YTick',      [0   5    10   15    20]); % Strike
set(gca,'YTickLabel', [95  97.5 100 102.5 105]);
alpha(0.5);
xlim([0 20]);
ylim([1 21]);
zlim([0 .3]);

% Surface plot for local volatility
hold on
surf(LV)
hold off

