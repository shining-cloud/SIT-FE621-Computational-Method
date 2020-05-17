function [pu pm pd] = probY(Vt,V0,Yt,dt,rho,sigma,kappa)

% Obtain the Y-tree probabilities for the Belieava-Nawalkha tree
% INPUTS
%   Vt = current value of V(t)
%   V0 = value of V(0) (Heston parameter)
%   dt = time increment
%   k  = ceiling(sqrt(V(t)/V(0))
%   rho = Heston CIR parameter
%   sigma = Heston CIR parameter
%   kappa = Heston CIR parameter
% OUTPUTS
%   Single probabilities pu,pm,pd

% Equation 11
if Vt > 0
    k = ceil(sqrt(Vt/V0));
else
    k = 1;
end

% Equations 9 and 6
muY = (rho/sigma*kappa - 0.5)*Vt;
sigmayt = sqrt(1-rho^2)*sqrt(Vt);
sigmay0 = sqrt(1-rho^2)*sqrt(V0);

% Calculate the Yu, Ym, Yd from Equation 13
I = (round(muY/k/sigmay0*sqrt(dt)));
Yu = Yt + (I+1)*k*sigmay0*sqrt(dt);
Ym = Yt + (I+0)*k*sigmay0*sqrt(dt);
Yd = Yt + (I-1)*k*sigmay0*sqrt(dt);

% Equation 17
eu = Yu - Yt - muY*dt;
em = Ym - Yt - muY*dt;
ed = Yd - Yt - muY*dt;

% Equation 16
pu = 0.5*(sigmayt^2*dt + em*ed) / k^2 / sigmay0^2 / dt;
pm =    -(sigmayt^2*dt + eu*ed) / k^2 / sigmay0^2 / dt;
pd = 0.5*(sigmayt^2*dt + eu*em) / k^2 / sigmay0^2 / dt;

% while (pd<0 || pm<0 || pu<0)
%     k = k+1;
%     I = (round(muY/k/sigmay0*sqrt(dt)));
%     Yu = Yt + (I+1)*k*sigmay0*sqrt(dt);
%     Ym = Yt + (I+0)*k*sigmay0*sqrt(dt);
%     Yd = Yt + (I-1)*k*sigmay0*sqrt(dt);
%     eu = Yu - Yt - muY*dt;
%     em = Ym - Yt - muY*dt;
%     ed = Yd - Yt - muY*dt;
%     pu = 0.5*(sigmayt^2*dt + em*ed) / k^2 / sigmay0^2 / dt;
%     pm =    -(sigmayt^2*dt + eu*ed) / k^2 / sigmay0^2 / dt;
%     pd = 0.5*(sigmayt^2*dt + eu*em) / k^2 / sigmay0^2 / dt;
% end
% 
    