function y = CZCharFun(S0,tau,t,params,K,rf,q,phi,psi,FunNum)

% Returns the characteristic function for the risk neutral probability P2.
% phi = integration variable
% psi = integration variable
% t  = maturity differential
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0     = initial variance
% Option features.
%    tau = maturity
%    K   = strike price
%    S   = spot price
%    rf  = risk free rate
%    q   = dividend yield

kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);

% Log of the stock price.
x = log(S0);

% Parameters "a" and "b"
a = kappa*theta;
b = kappa + lambda;

% "d" and "g" functions
d = sqrt((rho*sigma*i*phi - b)^2 + sigma^2*phi*(phi+i));
g = (b - rho*sigma*i*phi - sigma^2*i*psi + d) ...
  / (b - rho*sigma*i*phi - sigma^2*i*psi - d);

% The components of the affine characteristic function.
G = (1-g*exp(d*(tau-t)))/(1-g);

C = (rf-q)*i*phi*(tau-t) ...
  + a/sigma^2*((b - rho*sigma*i*phi + d)*(tau-t) - 2*log(G));
F = (1-exp(d*(tau-t)))/(1-g*exp(d*(tau-t)));
D = i*psi + (b - rho*sigma*i*phi - sigma^2*i*psi + d)/sigma^2 * F;

% The second characteristic function.
f2 = exp(C + D*v0 + i*phi*x);

if (FunNum == 2)
    y = f2;
else
    d = sqrt((rho*sigma*i*(phi-i) - b)^2 + sigma^2*(phi-i)*phi);
    g = (b - rho*sigma*i*(phi-i) - sigma^2*i*psi + d) ...
        / (b - rho*sigma*i*(phi-i) - sigma^2*i*psi - d);

    % The components of the affine characteristic function.
    G = (1-g*exp(d*(tau-t)))/(1-g);

    C = (rf-q)*i*(phi-i)*(tau-t) ...
        + a/sigma^2*((b - rho*sigma*i*(phi-i) + d)*(tau-t) - 2*log(G));
    F = (1-exp(d*(tau-t)))/(1-g*exp(d*(tau-t)));
    D = i*psi + (b - rho*sigma*i*(phi-i) - sigma^2*i*psi + d)/sigma^2 * F;

    % The second characteristic function with phi = phi-i
    F2 = exp(C + D*v0 + i*(phi-i)*x);

    % The first characteristic function
    y = 1/S0 * exp(-(rf-q)*(tau-t)) * F2;
end

