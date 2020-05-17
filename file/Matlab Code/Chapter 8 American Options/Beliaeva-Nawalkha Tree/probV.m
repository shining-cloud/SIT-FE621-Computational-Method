function [pu pm pd] = probV(X,X0,dt,kappa,theta,sigma)

% Obtain the V-tree probabilities for the Belieava-Nawalkha tree
% INPUTS
%   X = single transformed volatility
%   X0 = time zero transformed volatility
%   dt = time increment
%   k  = ceiling(sqrt(V(t)/V(0))
%   rho = Heston CIR parameter
%   sigma = Heston CIR parameter
%   kappa = Heston CIR parameter
% OUTPUTS
%  Matrix of probabilities pu,pm,pd

% Initial quantities
be = X0/sqrt(dt)/floor(X0/sqrt(1.5*dt));
bc = X0/sqrt(dt)/floor(X0/sqrt(1.5*dt)+1);
if abs(bc-sqrt(1.5)) < abs(be-sqrt(1.5))
	b = bc;
else
	b = be;
end

% Equations (22) and (30)
muX = 1/X*(0.5*kappa*(4*theta/sigma^2 - X^2) - 0.5);
J =  floor(muX*sqrt(dt)/b + 1/b^2);

if X > 0
    % Probabilities where X > 0 (Equation 28)
    pu = 1/2/b^2 - J/2 + 1/2/b*muX*sqrt(dt);
    pm = 1 - 1/b^2;
    pd = 1/2/b^2 + J/2 - 1/2/b*muX*sqrt(dt);
else
    % Probabilities where X = 0 (Equation 33)
    Xu = X + b*(J+1)*sqrt(dt);
    Vu = Xu^2*sigma^2/4;
    pu = kappa*theta*dt/Vu;
    pm = 0;
    pd = 1 - pu;
end

