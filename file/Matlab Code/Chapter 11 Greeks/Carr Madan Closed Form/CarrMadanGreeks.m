function y = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap,Greek,x,w)

% Perform the numerical integration
for j=1:length(x)
     u    = x(j);
     I(j) = w(j) * exp(-i*u*log(K)) * exp(-r*tau) ...
          * HestonCFGreek(u-(alpha+1)*i,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap,Greek) ...
          / (alpha^2 + alpha - u^2 + i*(2*alpha+1)*u);
     I(j) = real(I(j));
end

% Calculate the desired Greek
if strcmp(Greek,'Delta') || strcmp(Greek,'Gamma') || strcmp(Greek,'Rho') ||...
   strcmp(Greek,'Theta') || strcmp(Greek,'Volga') || strcmp(Greek,'Price');
    y = exp(-alpha*log(K))*sum(I)/pi;
elseif strcmp(Greek,'Vega1') || strcmp(Greek,'Vanna')
    y = exp(-alpha*log(K))*sum(I)/pi * 2*sqrt(v0);
end

