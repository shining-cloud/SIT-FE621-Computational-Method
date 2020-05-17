function y = LewisIntegrand(k,X,v0,tau,theta,kappa,sigma,rho)

% Integrand for the Lewis (2000) implementation of the Heston model, using
% the Fundamental Transform.

% k = Integration variable
% X = log(S/K) + (r-q)*tau;
% tau = moneyness
% v0, theta, kappa, sigma, rho = heston parameters

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
y = real(exp(-X*i*k)/(k^2 - i*k)*H);


