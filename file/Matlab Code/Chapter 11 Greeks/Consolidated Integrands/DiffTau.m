function [dC dD] = DiffTau(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap)

a = kappa*theta;
u = -0.5;
b = kappa + lambda;
d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);

if trap==1
    % "Little Heston Trap" formulation
    c = 1/g;
    D = (b - rho*sigma*i*phi - d)/sigma^2*((1-exp(-d*tau))/(1-c*exp(-d*tau)));
    G = (1 - c*exp(-d*tau))/(1-c);
    C = (r-q)*i*phi*tau + a/sigma^2*((b - rho*sigma*i*phi - d)*tau - 2*log(G));
elseif trap==0
    % Original Heston formulation.
    G = (1 - g*exp(d*tau))/(1-g);
    C = (r-q)*i*phi*tau + a/sigma^2*((b - rho*sigma*i*phi + d)*tau - 2*log(G));
    D = (b - rho*sigma*i*phi + d)/sigma^2*((1-exp(d*tau))/(1-g*exp(d*tau)));
end

% Derivatives of C and D w.r.t. maturity tau
dD = d*exp(d*tau)*(b-rho*sigma*phi*i+d)*(g-1)/sigma^2/(1-g*exp(d*tau))^2;
dC = (r-q)*phi*i + kappa*theta/sigma^2 ...
    * ((b-rho*sigma*phi*i+d) + 2*g*d*exp(d*tau)/(1-g*exp(d*tau)));

