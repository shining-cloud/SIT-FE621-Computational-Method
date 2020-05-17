function y = LordKahlFindAlpha(alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap)

PSI = HestonPsi(0,alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);

y = -alpha*log(K)  + (1/2)*log(PSI^2);

