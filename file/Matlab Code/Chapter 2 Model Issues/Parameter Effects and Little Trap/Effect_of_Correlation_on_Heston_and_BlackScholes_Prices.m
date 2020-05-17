% Difference between Heston and Black Scholes prices
% Effect of correlation

clc; clear;

% Black Scholes call
BSC = @(s,K,r,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-r*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));

kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.01;    % Heston parameter: mean reversion level.
v0 = 0.01;       % Heston parameter: initial variance.
lambda = 0;      % Heston parameter: risk.
sigma = 0.1;     % Heston parameter: volatility of variance.
K = 100;         % Spot price.
r = 0.0;         % Risk free rate.
q = 0.0;         % Dividend yield
tau = .5;        % Time to maturity.
PutCall = 'C';   % 'C'all or 'P'ut.
rhopos = 0.5;    % Heston parameter: correlation.
rhoneg = -0.5;
trap = 1;        % 1 = "Little Trap" form, 0 = Original Heston form
Lphi = 1e-15;    % Lower integration limit
dphi = 1e-2;     % Integration increment
Uphi = 100;      % Upper integration limit

% Black Scholes matched implied volatilities
volneg = 0.071037274323352 * sqrt(2);
vol0   = 0.070712338973920 * sqrt(2);
volpos = 0.070386797400082 * sqrt(2);

% Gauss Laguerre points
[x w] = GenerateGaussLaguerre(32);

% Generate the strikes.
S = [70:140];


%% Generate the Heston and the Black Scholes prices, and their difference
for i=1:length(S);
 	HCallpos(i) = HestonPriceTrapezoidal(S(i),K,tau,r,q,kappa,theta,sigma,lambda,v0,rhopos,trap,PutCall,dphi,Uphi,Lphi);
 	HCallneg(i) = HestonPriceTrapezoidal(S(i),K,tau,r,q,kappa,theta,sigma,lambda,v0,rhoneg,trap,PutCall,dphi,Uphi,Lphi);
	BSCallpos(i) = BSC(S(i),K,r,q,volpos,tau);
	BSCallneg(i) = BSC(S(i),K,r,q,volneg,tau);
	diffpos(i)   = HCallpos(i) - BSCallpos(i);	
	diffneg(i)   = HCallneg(i) - BSCallneg(i);	
end


%% Plot the price differences
plot(S,diffpos,'r+-',S,diffneg,'k--')
xlabel('Spot Price')
ylabel('Heston minus Black-Scholes')
legend('Rho = +0.5', 'Rho = -0.5');
axis([70 140 -0.15 .15])
