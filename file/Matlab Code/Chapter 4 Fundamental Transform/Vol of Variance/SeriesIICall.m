function [ivx y] = SeriesIICall(S,K,rf,q,T,v0,rho,theta,kappa,sigma)

% Series II call price (y) and (ivx)
% S = Spot Price
% K = strike price
% rf = risk free rate
% q = dividend yield
% T = maturity
% v0,rho,theta,kappa,sigma = Heston parameters


% Black Scholes call price and Black Scholes vega
% Note that v = variance, not volatility
BSC = @(S,K,rf,q,v,T) (S*exp(-q*T)*normcdf((log(S/K) + T*(rf - q + 0.5*v)) / sqrt(v*T)) - exp(-rf*T)*K*normcdf((log(S/K) + T*(rf - q + 0.5*v)) / sqrt(v*T) - sqrt(v*T)));

% The "J" integrals
J1 = J(rho,theta,kappa,T,v0,1);
J3 = J(rho,theta,kappa,T,v0,3);
J4 = J(rho,theta,kappa,T,v0,4);

% Time average of the deterministic volatility
v = theta + (v0-theta)*(1-exp(-kappa*T))/(kappa*T);

% X and Z required for the ratios of Black Scholes derivatives
X = log(S*exp((rf-q)*T)/K);
Z = v*T;

% The ratios of Black Scholes derivatives
R20 = R(2,0,X,Z,T);
R11 = R(1,1,X,Z,T);
R12 = R(1,2,X,Z,T);
R22 = R(2,2,X,Z,T);

% Series II volatility of volatility expansion for the implied volatility
iv = v ...
   + sigma/T*J1*R11 ...
   + sigma^2*(J3*R20/T^2 + J4*R12/T + 0.5/T^2*J1^2*(R22-R11^2*R20));

% Series II call price
y = BSC(S,K,rf,q,iv,T);

% Return the volatility
ivx = sqrt(iv);