function y = D(phi,kappa,theta,lambda,rho,sigma,tau,r,q,Trap)

a = kappa*theta;
u = -0.5;
b = kappa + lambda;
d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);

if Trap==1
	% "Little Heston Trap" formulation
    c = 1/g;
	y = (b - rho*sigma*i*phi - d)/sigma^2*((1-exp(-d*tau))/(1-c*exp(-d*tau)));
elseif Trap==0
	% Original Heston formulation.
	y = (b - rho*sigma*i*phi + d)/sigma^2*((1-exp(d*tau))/(1-g*exp(d*tau)));
end
