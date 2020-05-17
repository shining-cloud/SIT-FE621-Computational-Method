function [b0 b1] = findB(tau,params,K,r,q,V00,V10,b00,b10,xs,ws,xt,wt,Nt,Ntau,tol0,tol1,Ntol,a,b,c,d,DoubleType)

% Function to find b0 and b1 for the Chiarella and Ziogas American call
%   tau = maturity
%   params = Heston parameters
%   K   = strike price
%   rf  = risk free rate
%   q   = dividend yield
%   V00 = starting value, v0(n)
%   V10 = starting value, v1(n)
%   b00 = starting value, b0(n)
%   b10 = starting value, b1(n)
%   xs = abscissas for Gauss-Laguerre
%   ws = weights for Gauss-Laguerre
%   xt = abscissas for Gauss-Legendre
%   wt = weights for Gauss-Legendre
%   Nt = number of points for double trapezoidal rule for early exercise
%   Ntau = number of loops for the maturities
%   tol0 = tolerance on b0 for k-step 
%   tol1 = tolerance on b1 for k-step 
%   Ntol = tolerance for Newton's method
%   integration grid (t,x) on (a,b) x (c,d)
%   DoubleType = type of double integration for early exercise

% Conditional expected value of variance
EV = @(Vs,theta,kappa,t,s) (theta + (Vs - theta)*exp(-kappa*(t-s)));

% Initialize the quantities
v0 = zeros(Ntau,1);
v1 = zeros(Ntau,1);
b0 = zeros(Ntau,1);
b1 = zeros(Ntau,1);

% Starting values
v0(1) = V00;
v1(1) = V10;
b0(1) = b00;
b1(1) = b10;

% Parameters and time increment
kappa = params(1);
theta = params(2);
sigma = params(3);
v0t   = params(4);
rho   = params(5);
lambda = params(6);
dtau = tau/Ntau;

% First maturity step
fprintf('  n       mat        v0(n)      v1(n)      b0(n)      b1(n)  NumSteps \n')
fprintf('-------------------------------------------------------------------------\n')
fprintf('%4.0f %10.4f %10.4f %10.4f %10.4f %10.4f \n',1,dtau,v0(1),v1(1),b0(1),b1(1));

% Remaining maturity steps
for n=2:Ntau
    % Maturity
    mat = n*dtau;
    % Variances
    Evt(n) = EV(v0t,theta,kappa,tau,mat);
    v0(n)  = Evt(n) + sigma/kappa*sqrt(kappa*theta/2);
    v1(n)  = Evt(n) - sigma/kappa*sqrt(kappa*theta/2);
    % Starting values for Newton's method
    b0k_ = b0(n-1);
    b1k_ = b1(n-1);
    % Set the counter and the initial differences
    counter = 0;
    diff0 = 1.1*tol0;
    diff1 = 1.1*tol1;
    % Loop through until Newton's method returns bk and bk(-1) to tolerance
    while (diff0>tol0) && (diff1>tol1)
        counter = counter+1;
        % b1k = CZNewton(b1k_,v0(n),v1(n),tau,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k_,1,Ntol,a,b,c,d,DoubleType);
        % b0k = CZNewton(b0k_,v0(n),v1(n),tau,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k ,0,Ntol,a,b,d,c,DoubleType);
        b1k = CZNewton(b1k_,v0(n),v1(n),mat,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k_,1,Ntol,a,b,c,d,DoubleType);
        b0k = CZNewton(b0k_,v0(n),v1(n),mat,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k ,0,Ntol,a,b,d,c,DoubleType);
        b0(n) = b0k;
        b1(n) = b1k;
        diff0 = abs(b0k_ - b0k);
        diff1 = abs(b1k_ - b1k);
        b0k_ = b0k;
        b1k_ = b1k;
    end
    fprintf('%4.0f %10.4f %10.4f %10.4f %10.4f %10.4f %5.0f \n',n,mat,v0(n),v1(n),b0(n),b1(n),counter);
end
fprintf('-------------------------------------------------------------------------\n')

