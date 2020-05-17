% Plot of the Heston integrand
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
sigma1 = 0.175;    % Heston parameter: Volatility of variance
v01    = 0.01;     % Heston parameter: Initial variance
theta1 = 0.01;     % Heston parameter: Mean reversion level
tau1 = 1/52;       % Time to maturity
S1 = 7;            % Spot price

% Settings for the second integrand
sigma2 = 0.09;    % Heston parameter: Volatility of variance
v02    = 0.07;    % Heston parameter: Initial variance
theta2 = 0.07;    % Heston parameter: Mean reversion level
tau2 = 1;         % Time to maturity
S2 = 10;          % Spot price

% Integration grid
dphi = 0.1;       % Increment
Uphi = 100;       % Upper limit
Lphi = 0.00001;   % Lower limit
PHI = [Lphi:dphi:Uphi];

%% Calculate the integrands
for x=1:length(PHI)
	phi = PHI(x);
	Inte1(x) = HestonProb(phi,kappa,theta1,lambda,rho,sigma1,tau1,K,S1,r,q,v01,Pnum,Trap);
	Inte2(x) = HestonProb(phi,kappa,theta2,lambda,rho,sigma2,tau2,K,S2,r,q,v02,Pnum,Trap);
end

%% Plot the integrands
plot(PHI,Inte1,'k-',PHI,Inte2,'r-')
legend('First Integrand', 'Second Integrand')
xlabel('Integration range')
ylabel('Integrand')
ylim([-.1 .1])
