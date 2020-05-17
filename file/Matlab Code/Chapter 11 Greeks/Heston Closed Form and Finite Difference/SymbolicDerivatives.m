% Derivatives of C and D w.r.t. rho and sigma
clc; clear;

syms kappa theta sigma sigma rho i phi tau r q x v

% First characteristic function j=1
% Original Heston formulation.
u = 0.5;
b = kappa - rho*sigma;
d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
G = (1 - g*exp(d*tau))/(1-g);
C1 = (r-q)*i*phi*tau + kappa*theta/sigma^2*((b - rho*sigma*i*phi + d)*tau - 2*log(G));
D1 = (b - rho*sigma*i*phi + d)/sigma^2*((1-exp(d*tau))/(1-g*exp(d*tau)));

% Second characteristic function j=2
% Original Heston formulation.
u = -0.5;
b = kappa;
d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
G = (1 - g*exp(d*tau))/(1-g);
C2 = (r-q)*i*phi*tau + kappa*theta/sigma^2*((b - rho*sigma*i*phi + d)*tau - 2*log(G));
D2 = (b - rho*sigma*i*phi + d)/sigma^2*((1-exp(d*tau))/(1-g*exp(d*tau)));

% Derivatives with respect to rho
dC1dr = diff(C1,rho);
dD1dr = diff(D1,rho);
dC2dr = diff(C2,rho);
dD2dr = diff(D2,rho);

% Derivatives with respect to sigma
dC1ds = diff(C1,sigma);
dD1ds = diff(D1,sigma);
dC2ds = diff(C2,sigma);
dD2ds = diff(D2,sigma);

% Derivatives with respect to kappa
dC1dk = diff(C1,kappa)
dD1dk = diff(D1,kappa)
dC2dk = diff(C2,kappa)
dD2dk = diff(D2,kappa)


