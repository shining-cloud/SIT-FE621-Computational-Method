function y = LewisIntegrandGreek(k,S,r,q,X,v0,tau,theta,kappa,sigma,rho,Greek)
% Integrand for the Lewis (2000) implementation of the Heston model, using
% the Fundamental Transform.

% Define the updated "tilde" parameters
kappa = 2*(kappa + i*k*rho*sigma)/sigma^2;
theta = kappa*theta/(kappa + i*k*rho*sigma);
t = tau*sigma^2/2;
c = (k^2-i*k)/sigma^2;

% Compute the required quantities
d = sqrt(kappa^2 + 4*c);
alpha = (-kappa + d)/2;
beta  = (-kappa - d)/2;
g = beta/alpha;

% The functions B(t) and A(t)
D = (kappa+d)*(1-exp(d*t)) / (1 - g*exp(d*t))/2;
C = kappa*theta*((kappa+d)*t/2 - log((1-g*exp(d*t))/(1-g)) );

% The fundamental transform, H(k,S,t)
H = exp(C + D*v0);

% Return the real part of the integrand.
if strcmp(Greek,'Price')
    y = real(exp(-X*i*k-r*tau)/(k^2-i*k)*H);
elseif strcmp(Greek,'Delta')
    y = real(exp(-X*i*k-r*tau)/(k^2-i*k)*H*(-i*k/S));
elseif strcmp(Greek,'Gamma')
    y = -real(exp(-X*i*k-r*tau)*H/S^2);
elseif strcmp(Greek,'Vega1')
    y = real(exp(-X*i*k-r*tau)/(k^2-i*k)*H*D)*2*sqrt(v0);
elseif strcmp(Greek,'Rho')
    y = real(exp(-X*i*k-r*tau)*(-i*k*tau-tau)/(k^2-i*k)*H);
elseif strcmp(Greek,'Theta')
    dC = kappa*theta*((kappa+d)/2 + g*d*exp(d*t)/(1-g*exp(d*t))) * sigma^2/2;
    dD = (kappa+d)/2*exp(d*t)*d*(g-1)/(1-g*exp(d*t))^2 * sigma^2/2;
    y =  real((-i*k*(r-q)-r + (dC+dD*v0))*H*exp(-X*i*k-r*tau)/(k^2-i*k));
elseif strcmp(Greek,'Volga')
    I = exp(-X*i*k-r*tau)/(k^2-i*k)*H*D*D*sqrt(v0) ...
      + exp(-X*i*k-r*tau)/(k^2-i*k)*H*D/2/sqrt(v0);
    y = real(I)*4*sqrt(v0); 
elseif strcmp(Greek,'Vanna')
    y = real(exp(-X*i*k-r*tau)/(k^2-i*k)*H*D*(-i*k/S))*2*sqrt(v0);
end

