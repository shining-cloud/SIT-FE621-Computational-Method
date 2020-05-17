% Effect of volatility of variance (sigma) on Heston prices relative to
% Black Scholes

clc; clear;

% Black Scholes call
BSC = @(s,K,r,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-r*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));

v0 = 0.01;       % Heston parameter: initial variance.
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.01;    % Heston parameter: mean reversion level.
lambda = 0;      % Heston parameter: risk.
K = 100;         % Spot price.
r = 0.0;         % Risk free rate.
q = 0.0;         % Dividend Yield
tau = .5;        % Time to maturity.
PutCall = 'C';   % 'C'all or 'P'ut.
trap = 1;        % 1 = "Little Trap" form, 0 = Original Heston form
Lphi = 1e-15;    % Trapezoidal rule Lower integration limit
dphi = 1e-2;     % Trapezoidal rule Integration increment
Uphi = 100;      % Trapezoidal rule Upper integration limit

rho = 0;         % Heston parameter rho

% Black Scholes matched implied volatilities
volneg = 0.071037274323352 * sqrt(2);
vol0   = 0.070712338973920 * sqrt(2);
volpos = 0.070386797400082 * sqrt(2);

% Take the matched volatilty corresponding to rho = 0
vol = vol0;

% Two values of sigma
sigmalo = 0.1;
sigmahi = 0.2;

% Generate the strikes.
S = [70:0.25:140];

% Generate Gauss Laguerre abscissas and weights
[x w] = GenerateGaussLaguerre(32);

%% Generate the Heston prices, and the implied volatilities from those prices.
for i=1:length(S);
    HCallhi(i) = HestonPriceTrapezoidal(S(i),K,tau,r,q,kappa,theta,sigmahi,lambda,v0,rho,trap,PutCall,dphi,Uphi,Lphi);
    HCalllo(i) = HestonPriceTrapezoidal(S(i),K,tau,r,q,kappa,theta,sigmalo,lambda,v0,rho,trap,PutCall,dphi,Uphi,Lphi);
	BSCall(i)  = BSC(S(i),K,r,q,vol,tau);
	diffhi(i)  = HCallhi(i) - BSCall(i);	
	difflo(i)  = HCalllo(i) - BSCall(i);	
end


%% Plot the implied volatility skew.
plot(S,difflo,'k-',S,diffhi,'r-')
xlabel('Spot Price')
ylabel('Heston minus Black-Scholes')
legend('Sigma = 0.1', 'Sigma = 0.2');

