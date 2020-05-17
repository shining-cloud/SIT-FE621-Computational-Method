function C = Ct(phi,kappa,theta,sigma,rho,rf,q,T,Pnum,C0,D0)

% Calculates the C(t) coefficient for the M-N time dependent Heston model
% phi = integration variable
% kappa,theta,sigma,rho = Heston parameters
% rf = risk free rate
% q = dividend yield
% T = maturity
% Pnum = characteristic function number (1 or 2)
% C0 = previous value of C(t-1)
% D0 = previous value of D(t-1)


lambda = 0;
if Pnum==1
	u = 0.5;
	b = kappa + lambda - rho*sigma;
else
	u = -0.5;
	b = kappa + lambda;
end

d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d - D0*sigma^2) / (b - rho*sigma*i*phi - d - D0*sigma^2);
G = (1 - g*exp(d*T))/(1-g);
C = (rf-q)*i*phi*T + kappa*theta/sigma^2*((b - rho*sigma*i*phi + d)*T - 2*log(G)) + C0;
