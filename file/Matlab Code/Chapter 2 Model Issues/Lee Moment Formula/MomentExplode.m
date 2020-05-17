function y = MomentExplode(w,lambda,sigma,kappa,rho)

% Andersen and Piterbarg moment explosion
% INPUTS
%   w = moment number
%   lambda = Heston drift (set to 1)
%   sigma  = Heston vol of variance
%   kappa  = Hestn mean reversion speed
%   rho    = Heston correlation
% OUTPUT
%   T = Time of moment explosion
%   T < inf --> moments explode
%   T = inf --> moments don't explode

k = lambda^2*w*(w-1)/2;
b = 2*k/sigma^2;
a = 2*(rho*sigma*lambda*w - kappa)/sigma^2;
D = a^2 - 4*b;
if D>=0 & a<0
    T = inf;
elseif D>=0 & a>0
    g = sqrt(D)/2;
    T = log((a/2+g)/(a/2-g))/g/sigma^2;
elseif D<0
    beta = sqrt(-D)/2;
    if a<0
        PI = pi;
    else
        PI = 0;
    end
    T = 2*(PI + atan(2*beta/a))/beta/sigma^2;
end

y = T;

