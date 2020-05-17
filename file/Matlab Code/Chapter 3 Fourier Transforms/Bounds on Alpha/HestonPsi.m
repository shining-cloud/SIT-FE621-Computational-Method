function y = HestonPsi(v,alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap)

% "psi" function for Carr-Madan (1999) call price formulation
% v = integration variable
% alpha = damping factor
% kappa, theta, lambda, rho, sigma, v0 = Heston parameters
% tau = maturity
% K = strike price
% S = spot price
% r = risk free rate
% trap = 0 Little Trap, 1 Heston

% The Heston characteristic function
CF = HestonCF(v-(alpha+1)*i, kappa, theta, lambda, rho, sigma, tau, K, S, r, v0, trap);

% The Carr Madan psi function
y = exp(-r*tau) * exp(-i*v*log(K)) * CF ...
  / (alpha^2 + alpha - v^2 + i*v*(2*alpha+1));

y = real(y);


