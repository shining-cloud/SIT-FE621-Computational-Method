function [Price LoDomain HiDomain Npoints] = HestonPriceGaussLegendreMD(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLe,wGLe,A,tol)

% Multi-domain integration of Zhu (2010).
% Uses Gauss-Legendre abscissas and weights on each sub-domain
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T = Time to maturity.
%   r = Risk free rate.
%   kappa  = Heston parameter: mean reversion speed.
%   theta  = Heston parameter: mean reversion level.
%   sigma  = Heston parameter: volatility of vol
%   lambda = Heston parameter: risk.
%   v0     = Heston parameter: initial variance.
%   rho    = Heston parameter: correlation
%   trap:  1 = "Little Trap" formulation
%          0 = Original Heston formulation
%   xGLe = Gauss Legendre abscissas
%   wGLe = Gauss Legendre weights
%   A   = Array of points (0,a1], (a1,a2], (a2,a3] for the subdomain 
%   tol = Tolerance for stopping the subdomains
% OUTPUT -------------------------------------------------------
%   Price = The Heston call or put price
%   Domain = Resulting domain of integration
%   Npoints = Number of integration points

for j=2:length(A);
    for k=1:length(xGLe);
        % Lower and upper limits of the subdomain
        a = A(j-1);
        b = A(j);
        X = (a+b)/2 + (b-a)/2*xGLe(k);
        % The integrals
        int1(j,k) = wGLe(k)*HestonProb(X,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap)*(b-a)/2;
        int2(j,k) = wGLe(k)*HestonProb(X,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap)*(b-a)/2;
        % Sum the integrals over each subdomain
        sum1(j) = sum(int1(j,:));
        sum2(j) = sum(int2(j,:));
    end
    % Stopping criterion
    if abs(sum1(j))<tol && abs(sum2(j))<tol
        break;
    end
end

% Define P1 and P2
P1 = 1/2 + 1/pi*sum(sum1);
P2 = 1/2 + 1/pi*sum(sum2);

% The call price
HestonC = S*exp(-q*T)*P1 - K*exp(-r*T)*P2;

% The put price by put-call parity
HestonP = HestonC - S*exp(-q*T) + K*exp(-r*T);

% Output the option price
if strcmp(PutCall,'C')
	Price = HestonC;
else
	Price = HestonP;
end

% Output the integration domain
LoDomain = A(1);
HiDomain = A(j);
Domain = [A(1) A(j)];

% Output the number of points
Npoints = length(A(1:j))*length(xGLe);

