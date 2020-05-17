function y = CarrMadanIntegrand(u,alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap,PutCall);

% Returns the Carr-Madan integrand,
% based on the Heston characteristic function, f2 (the second CF)

if strcmp(PutCall,'C')
	I = exp(-i*u*log(K)) * exp(-r*tau) ...
		* HestonCF(u-(alpha+1)*i,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap) ...
		/ (alpha^2 + alpha - u^2 + i*(2*alpha+1)*u);
else
	I = exp(-i*u*log(K)) * exp(-r*tau) ...
		* HestonCF(u-(-alpha+1)*i,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap) ...
		/ (alpha^2 - alpha - u^2 + i*(-2*alpha+1)*u);
end
y = real(I);
