function y = CarrMadanDampedIntegrandOTM(v,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap,alpha)

% Returns the Carr-Madan integrand for OTM options
% based on the Heston characteristic function, f2 (the second CF)
% Calculate z(v-ia)
u  =  v - i*alpha;
phi = HestonCF(u-i,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap);
z1  = exp(-r*tau)*(S^(i*u+1)/(1+i*u) - exp(r*tau)*S^(i*u+1)/(i*u) - phi/(u^2-i*u));

% Calculate z(v+ia)
clear u phi
u  =  v + i*alpha;
phi = HestonCF(u-i,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap);
z2  = exp(-r*tau)*(S^(i*u+1)/(1+i*u) - exp(r*tau)*S^(i*u+1)/(i*u) - phi/(u^2-i*u));

% Calculate the Fourier transform of y 
y = exp(-i*u*log(K)) * (z1 - z2)/2;

% Return the real part only
y = real(y);
