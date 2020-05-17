% Plot of the Heston integrand along with Gauss Laguerre integration points
clc; clear;

K = 10;         % Strike Price
r = 0.0;        % Risk free rate
q = 0.0;        % Dividend yield
Pnum = 1;       % Characteristic function (j=1,2)
Trap = 0;       % Heston formulation for characteristic function

rho = -0.9;     % Heston parameter: correlation
lambda = 0;     % Heston parameter: risk preference
kappa = 10;     % Heston parameter: mean reversion speed

% Settings for the first integrand
sigma = 0.175;    % Heston parameter: Volatility of variance
v0    = 0.01;     % Heston parameter: Initial variance
theta = 0.01;     % Heston parameter: Mean reversion level
tau = 1/52;       % Time to maturity
S = 7;            % Spot price

% Integration grid
dphi = 0.1;       % Increment
Uphi = 120;       % Upper limit
Lphi = 1e-6;      % Lower limit
PHI = [Lphi:dphi:Uphi];

%% Calculate the integrand along the continuous grid
for k=1:length(PHI)
	phi = PHI(k);
	Int1(k) = HestonProb(phi,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Pnum,Trap);
end

%% Calculate the integrand along the Gauss-Laguerre points only
[x w] = GenerateGaussLaguerre(32);
for k=1:length(x)
    GLa(k) = HestonProb(x(k),kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Pnum,Trap);
end

%% Plot the integrands
plot(PHI,Int1,'k-',x,GLa,'bo')
legend('Integrand','Gauss-Laguerre points')
xlabel('Integration range')
ylabel('Integrand')
%ylim([-.1 .1])

