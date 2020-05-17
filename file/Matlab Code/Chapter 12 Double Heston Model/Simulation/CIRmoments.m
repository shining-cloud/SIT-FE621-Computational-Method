function [u1 u2] = CIRmoments(param,Vs,dt)

kappa = param(1);
theta = param(2);
sigma = param(3);

% E[vt | vs];
u1 = theta + (Vs - theta)*exp(-kappa*dt);

% Var[vt | vs]
s2 = Vs*sigma^2*exp(-kappa*dt)/kappa*(1-exp(-kappa*dt)) ...
  + theta*sigma^2/2/kappa*(1-exp(-kappa*dt))^2;

% E[vt^2 | vs]
u2 = s2 + u1^2;


