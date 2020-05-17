function y = A(phi1,phi2,param,tau,r,q,trap)

% Returns the Heston bivariate characteristic function coefficient "B"
%   phi1,phi2 = integration variables
%   param = Heston parameters:
% Option features.
%   tau = maturity
%   r = risk free rate
%   q = dividend yield
%   trap  1 = Little Trap formulation
%         0 = Original Heston formulation

kappa = param(1);
theta = param(2);
sigma = param(3);
rho   = param(4);
lambda = 0;

% Parameter "a" is the same for P1 and P2.
a = kappa*theta;

% Parameters "u" and "b" for P2.
u = -0.5;
b = kappa + lambda;
d = sqrt((rho*sigma*i*phi1 - b)^2 - sigma^2*(2*u*i*phi1 - phi1^2));

% Ensure denominator for g remains > 0
num = (b - rho*sigma*i*phi1 + d - sigma^2*i*phi2);
den = (b - rho*sigma*i*phi1 - d - sigma^2*i*phi2);
if den == 0
    den = 0.01;
end
g = num/den;
c = 1/g;

if trap==1
	% Little Trap formulation in Kahl (2008)
	G = (c*exp(-d*tau)-1)/(c-1);
	AA = (r-q)*i*phi1*tau + a/sigma^2*((b - rho*sigma*i*phi1 - d)*tau - 2*log(G));
elseif trap==0
	% Original Heston formulation.
	G = (1-g*exp(d*tau))/(1-g);
	AA = (r-q)*i*phi1*tau + a/sigma^2*((b - rho*sigma*i*phi1 + d)*tau - 2*log(G));
end

y = AA;
