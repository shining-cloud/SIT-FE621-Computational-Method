function y = B(phi,param,tau,trap,j)

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

if trap==1
	d1 = sqrt((kappa1-rho1*sigma1*i*phi)^2 + sigma1^2*phi*(phi+i));
	d2 = sqrt((kappa2-rho2*sigma2*i*phi)^2 + sigma2^2*phi*(phi+i));
	g1 = (kappa1-rho1*sigma1*phi*i-d1) / (kappa1-rho1*sigma1*phi*i+d1);
	g2 = (kappa2-rho2*sigma2*phi*i-d2) / (kappa2-rho2*sigma2*phi*i+d2);
	B1 = (kappa1-rho1*sigma1*phi*i-d1)*(1-exp(-d1*tau)) / sigma1^2 / (1-g1*exp(-d1*tau));
	B2 = (kappa2-rho2*sigma2*phi*i-d2)*(1-exp(-d2*tau)) / sigma2^2 / (1-g2*exp(-d2*tau));
else
	d1 = sqrt((kappa1-rho1*sigma1*phi*i)^2 + sigma1^2*(phi*i+phi^2));
	d2 = sqrt((kappa2-rho2*sigma2*phi*i)^2 + sigma2^2*(phi*i+phi^2));
	g1 = (kappa1-rho1*sigma1*phi*i+d1)/(kappa1-rho1*sigma1*phi*i-d1);
	g2 = (kappa2-rho2*sigma2*phi*i+d2)/(kappa2-rho2*sigma2*phi*i-d2);
	B1 = (kappa1-rho1*sigma1*phi*i+d1)*(1-exp(d1*tau))/sigma1^2/(1-g1*exp(d1*tau));
	B2 = (kappa2-rho2*sigma2*phi*i+d2)*(1-exp(d2*tau))/sigma2^2/(1-g2*exp(d2*tau));
end

if j==1
    y = B1;
else
    y = B2;
end
