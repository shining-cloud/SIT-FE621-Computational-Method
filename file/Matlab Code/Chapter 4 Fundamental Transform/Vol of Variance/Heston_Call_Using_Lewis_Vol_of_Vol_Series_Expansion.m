% Volatility of volatility expansion for the Heston call price
% from Alan Lewis "Option Valuation Under Stochastic Volatility" (2001)

clc; clear;

%% Parameter and other settings
% Spot price, strike price, risk free rate, dividend yield, and maturity
S   = 100;
K   = 105;
rf  = 0.05;
q   = 0.01;
T   = 0.25;

% Heston model parameters.
kappa = 10;      % Volatility reversion speed
theta = 0.07;    % Volatility reversion level
sigma = 0.9;     % Volatility of variance
rho   = 0.9;     % Correlation
v0    = 0.06;    % Initial variance
lambda = 0;      % Risk
trap = 1;        % Little trap formulation
					 
% Specify Newton Coates method
%   1 = mid-point rule, 2 = trapezoidal rule
%   3 = Simpson's rule, 4 = Simpson's 3/8 rule
method = 4;


%% Find the Heston prices using exact method and series expansions
% Find the exact Heston call price
a = 1e-10;
b = 100;
N = 10000;
ExactPrice = HestonPriceNewtonCoates('C',S,K,T,rf,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);

% Find the Series I Heston call price
SeriesIPrice = SeriesICall(S,K,rf,q,T,v0,rho,theta,kappa,sigma);

% Find the Series II implied volatility and call price
[IVx SeriesIIPrice] = SeriesIICall(S,K,rf,q,T,v0,rho,theta,kappa,sigma);

% Time average of the deterministic variance
v = theta + (v0-theta)*(1-exp(-kappa*T))/(kappa*T);

% Time average of the deterministic volatility
IV = sqrt(v);

% Black Scholes call price
% Note that v = variance, not volatility
BSC = @(S,K,rf,q,v,T) (S*exp(-q*T)*normcdf((log(S/K) + T*(rf - q + 0.5*v)) / sqrt(v*T)) - exp(-rf*T)*K*normcdf((log(S/K) + T*(rf - q + 0.5*v)) / sqrt(v*T) - sqrt(v*T)));
		 
% Use time average variance for the Black Scholes call price
BSPrice = BSC(S,K,rf,q,v,T);

%% Output the prices and IVs
clc;
fprintf('Method               Price \n');
fprintf('---------------------------\n');
fprintf('Black Scholes    %10.4f\n',BSPrice);
fprintf('Exact            %10.4f\n',ExactPrice);
fprintf('Series I         %10.4f\n',SeriesIPrice);
fprintf('Series II        %10.4f\n',SeriesIIPrice);
fprintf('---------------------------\n');
fprintf('Time average of deterministic vol = %6.4f\n',IV);
fprintf('Series II implied volatility      = %6.4f\n',IVx);
fprintf('Series II implied variance        = %6.4f\n',IVx^2);


