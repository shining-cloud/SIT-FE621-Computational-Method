function y = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,Greek)

lambda = 0;

% Compute the integral.
% The integration variable is complex k = kr + (1/2)*i
for j=1:length(x);
	u = x(j);
    int(j) = w(j)*LewisIntegrand311Greeks(u,kappa,theta,lambda,rho,sigma,T,S,K,r,q,v0,trap,Greek);
end
Integral = sum(int);

% Equation (3.11) in Lewis (2011)
if strcmp(Greek,'Price')
    y = S*exp(-q*T) - sqrt(K)/pi*Integral;
elseif strcmp(Greek,'Delta')
    y = exp(-q*T) - sqrt(K)/pi*Integral;
elseif strcmp(Greek,'Gamma') || strcmp(Greek,'Rho')
    y = -sqrt(K)/pi*Integral;
elseif strcmp(Greek,'Theta')
    y = q*S*exp(-q*T) + sqrt(K)/pi*Integral;
elseif strcmp(Greek,'Vega1') || strcmp(Greek,'Vanna')
    y = - sqrt(K)/pi*Integral*2*sqrt(v0);
elseif strcmp(Greek,'Volga')
    y = - sqrt(K)/pi*Integral;
end

