function y = DoubleHestonCF(phi,param,tau,S,rf,q,trap)

% Returns the integrand for the risk neutral probabilities P1 and P2.
% phi = integration variable
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v      = initial variance
% Option features.
%    S  = spot price
%    rf = risk free rate
%    trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

% First set of parameters
kappa1 = param(1);
theta1 = param(2);
sigma1 = param(3);
v01    = param(4);
rho1   = param(5);

% Second set of parameters
kappa2 = param(6);
theta2 = param(7);
sigma2 = param(8);
v02    = param(9);
rho2   = param(10);

x0 = log(S);

if trap==1
	d1 = sqrt((kappa1-rho1*sigma1*i*phi)^2 + sigma1^2*phi*(phi+i));
	d2 = sqrt((kappa2-rho2*sigma2*i*phi)^2 + sigma2^2*phi*(phi+i));
	G1 = (kappa1-rho1*sigma1*phi*i-d1) / (kappa1-rho1*sigma1*phi*i+d1);
	G2 = (kappa2-rho2*sigma2*phi*i-d2) / (kappa2-rho2*sigma2*phi*i+d2);
	B1 = (kappa1-rho1*sigma1*phi*i-d1)*(1-exp(-d1*tau)) / sigma1^2 / (1-G1*exp(-d1*tau));
	B2 = (kappa2-rho2*sigma2*phi*i-d2)*(1-exp(-d2*tau)) / sigma2^2 / (1-G2*exp(-d2*tau));
	X1 = (1-G1*exp(-d1*tau))/(1-G1);
	X2 = (1-G2*exp(-d2*tau))/(1-G2);
	A  = (rf-q)*phi*i*tau ...
		+ kappa1*theta1/sigma1^2*((kappa1-rho1*sigma1*phi*i-d1)*tau - 2*log(X1)) ...
		+ kappa2*theta2/sigma2^2*((kappa2-rho2*sigma2*phi*i-d2)*tau - 2*log(X2)) ;
else
	d1 = sqrt((kappa1-rho1*sigma1*phi*i)^2 + sigma1^2*(phi*i+phi^2));
	d2 = sqrt((kappa2-rho2*sigma2*phi*i)^2 + sigma2^2*(phi*i+phi^2));
	g1 = (kappa1-rho1*sigma1*phi*i+d1)/(kappa1-rho1*sigma1*phi*i-d1);
	g2 = (kappa2-rho2*sigma2*phi*i+d2)/(kappa2-rho2*sigma2*phi*i-d2);
	B1 = (kappa1-rho1*sigma1*phi*i+d1)*(1-exp(d1*tau))/sigma1^2/(1-g1*exp(d1*tau));
	B2 = (kappa2-rho2*sigma2*phi*i+d2)*(1-exp(d2*tau))/sigma2^2/(1-g2*exp(d2*tau));
	X1 = (1-g1*exp(d1*tau))/(1-g1);
	X2 = (1-g2*exp(d2*tau))/(1-g2);
	A = (rf-q)*phi*i*tau ...
		+ kappa1*theta1/sigma1^2*((kappa1-rho1*sigma1*phi*i+d1)*tau - 2*log(X1)) ...
		+ kappa2*theta2/sigma2^2*((kappa2-rho2*sigma2*phi*i+d2)*tau - 2*log(X2));
end

% The characteristic function.
y = exp(A + i*phi*x0 + B1*v01 + B2*v02);

