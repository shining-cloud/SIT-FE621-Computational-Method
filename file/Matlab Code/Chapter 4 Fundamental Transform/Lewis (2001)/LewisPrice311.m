function y = LewisPrice311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w)

% Implements the Heston model using the characteristic function described in
% the article by Alan Lewis (2001) "A Simple Option Formula for General 
% Jump Diffusion and Other Exponential Levy Processes"

% S = spot price
% K = strike price
% r = risk free rate
% q = dividend yield
% T = maturity
% theta, kappa, sigma, rho, v0 = Heston parameters
% trap = Little Trap (0) or Heston (1)
% x, w = abscissas and weights

lambda = 0;

% Compute the integral.
% The integration variable is complex k = kr - (1/2)*i
for j=1:length(x);
	u = x(j) - (1/2)*i;
    int(j) = w(j)*LewisIntegrand311(u,kappa,theta,lambda,rho,sigma,T,S,K,r,q,v0,trap);
end
Integral = sum(int);

% Equation (3.11) in Lewis (2011)
y = S*exp(-q*T) - (1/pi)*sqrt(K*S)*exp(-(r+q)*T/2)*Integral;
