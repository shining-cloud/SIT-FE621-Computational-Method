function y = dPjdT(Pnum,phi,kappa,theta,lambda,rho,sigma,T,K,S,rf,v0,Trap)

% Returns the integrand for the partial derivatives of P1 and P2 
% with respect to maturity, namely dP1/dT and dP2/dT.
% phi = integration variable
% Pnum = 1 or 2 (for the probabilities)
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0     = initial variance
% Option features.
%    T  = maturity
%    K  = strike price
%    S  = spot price
%    rf = risk free rate
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

% Log of the stock price.
x = log(S);

% Parameters "u" and "b" are different for P1 and P2.
if Pnum==1
	u = 0.5;
	b = kappa + lambda - rho*sigma;
else
	u = -0.5;
	b = kappa + lambda;
end

d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
edT = exp(d*T);

if Trap==1
	% Little Trap formulation
	c = 1/g;
	edT = exp(-d*T);	
	dDdT = (b - rho*sigma*i*phi - d)/sigma^2 ...
		* (d*edT*(1-c*edT) - (1-edT)*c*d*edT) ...
		/ (1 - c*edT)^2;
	dCdT = rf*i*phi ...
		+ kappa*theta/sigma^2 ...
		* ((b - rho*sigma*i*phi - d) - 2*c*d*edT/(1 - c*edT));
	G = (1 - c*edT)/(1-c);
	D = (b - rho*sigma*i*phi - d)/sigma^2*((1-edT)/(1-c*edT));
	C = rf*i*phi*T + kappa*theta/sigma^2*((b - rho*sigma*i*phi - d)*T - 2*log(G));
else
	% Original Heston formulation.
	edT = exp(d*T);
	dDdT = (b - rho*sigma*i*phi + d)/sigma^2 ...
		* (d*edT*(g*edT-1) + (1-edT)*g*d*edT) ...
		/ (1 - g*edT)^2;
	dCdT = rf*i*phi ...
		+ kappa*theta/sigma^2 ...
		* ((b - rho*sigma*i*phi + d) + 2*g*d*edT/(1 - g*edT));
	G = (1 - g*edT)/(1-g);
	C = rf*i*phi*T + kappa*theta/sigma^2*((b - rho*sigma*i*phi + d)*T - 2*log(G));
	D = (b - rho*sigma*i*phi + d)/sigma^2*((1-edT)/(1-g*edT));
end
	
dfdT = exp(C + D*v0 + i*phi*x)...
  	    * (dCdT + dDdT*v0);

y = real(K^(-i*phi) * dfdT / i / phi);
