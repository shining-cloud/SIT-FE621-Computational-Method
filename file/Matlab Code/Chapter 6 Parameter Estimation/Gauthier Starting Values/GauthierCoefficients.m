function [A B C D] = GauthierCoefficients(kappa,theta,v0,S,K,T,rf,q)

% Returns the Gauthier and Rivaille coefficienta
% kappa,theta,v0 = Heston parameters
% S = Spot price
% K = Strike price
% T = Maturity
% rf = Risk free rate
% q = Dividend yield

x = log(S);

kT = kappa*T;

m0 = exp(-kT)*(exp(kT)-1)/kappa;
m1 = T - m0;
varT = m0*v0 + m1*theta;

p0 = exp(-kT)*(exp(kT) - kT - 1)/kappa^2;
p1 = exp(-kT)*(exp(kT)*(kT-2) + kT + 2)/kappa^2;

q0 = exp(-kT)*(2*exp(kT) - kT*(kT+2) - 2)/2/kappa^3;
q1 = exp(-kT)*(2*exp(kT)*(kT-3) + kT*(kT+4) + 6)/2/kappa^3;

r0 = exp(-2*kT)*(2*exp(2*kT) - 4*exp(kT)*kT - 2)/4/kappa^3;
r1 = exp(-2*kT)*(exp(2*kT)*(2*kT-5) + 4*exp(kT)*(kT+1) + 1)/4/kappa^3;

% Generate the Black-Scholes derivatives
[P11 P21 P02 P22] = BlackScholesDerivatives(kappa,theta,v0,S,K,T,rf,q);
 
y = varT;

% Black Scholes Put Price
g = y^(-1/2)*(-x + log(K) - (rf-q)*T) - (1/2)*sqrt(y);
f = y^(-1/2)*(-x + log(K) - (rf-q)*T) + (1/2)*sqrt(y);
BSPut = K*exp(-rf*T)*normcdf(f) - S*exp(-q*T)*normcdf(g);

% Return the coefficients
A = BSPut;
B = (v0*r0 + theta*r1)*P02;
C = (v0*p0 + theta*p1)*P11;
D = (v0*q0 + theta*q1)*P21 + 0.5*(v0*p0 + theta*p1)^2*P22;

