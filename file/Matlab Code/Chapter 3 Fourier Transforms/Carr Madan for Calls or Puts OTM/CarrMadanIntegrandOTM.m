function z = CarrMadanIntegrandOTM(u,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap)

% Returns the Carr-Madan integrand for OTM options
% based on the Heston characteristic function, f2 (the second CF)

% The Heston characteristic function
phi = HestonCF(u-i,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap);

% Calculate the integrand
I = exp(-i*u*log(K)) * exp(-r*tau) ...
  * (S^(i*u+1)/(1+i*u) - exp(r*tau)*S^(i*u+1)/(i*u) - phi/(u^2-i*u));

% Return the real part
z = real(I);

