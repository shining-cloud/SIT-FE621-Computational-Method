function [dA dB1 dB2] = DiffTau(phi,param,tau,rf,q)

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

d1 = sqrt((kappa1-rho1*sigma1*phi*i)^2 + sigma1^2*(phi*i+phi^2));
d2 = sqrt((kappa2-rho2*sigma2*phi*i)^2 + sigma2^2*(phi*i+phi^2));
g1 = (kappa1-rho1*sigma1*phi*i+d1)/(kappa1-rho1*sigma1*phi*i-d1);
g2 = (kappa2-rho2*sigma2*phi*i+d2)/(kappa2-rho2*sigma2*phi*i-d2);

% The derivatives
dC1 = kappa1*theta1/sigma1^2 * ((kappa1-rho1*sigma1*phi*i+d1) + 2*g1*d1*exp(d1*tau)/(1-g1*exp(d1*tau)));
dC2 = kappa2*theta2/sigma2^2 * ((kappa2-rho2*sigma2*phi*i+d2) + 2*g2*d2*exp(d2*tau)/(1-g2*exp(d2*tau)));

dA = (rf-q)*phi*i + dC1 + dC2;
dB1 = d1*exp(d1*tau)*(kappa1-rho1*sigma1*phi*i+d1)*(g1-1)/sigma1^2/(1-g1*exp(d1*tau))^2;
dB2 = d2*exp(d2*tau)*(kappa2-rho2*sigma2*phi*i+d2)*(g2-1)/sigma2^2/(1-g2*exp(d2*tau))^2;

