function [y Int] = DoubleTrapezoidal(params,S0,K,tau,rf,q,b0,b1,X,T,funNum)

% Trapezoidal rule for double integrals
% params = Heston parameters
% S0  = spot price
% K   = strike price
% tau = maturity
% rf  = risk free interest rate
% q   = dividend yield
% b0  = Chiarella and Ziogas parameter
% b1  = Chiarella and Ziogas parameter
% T = vector of grid points for maturity direction
% X = vector of grid points for stock price direction
% funNum = 1 for f1 c.f., and 2 for f2 c.f.

if funNum == 1;
    rq = q;
elseif funNum == 2;
    rq = rf;
end

y = 0;
for t=2:length(T)
    a = T(t-1);
    b = T(t);
    for x=2:length(X)
        c = X(x-1);
        d = X(x);
        g1 = exp(rq*a) * real(exp(-b0*i*c) * CZCharFun(S0,tau,a,params,K,rf,q,c,-b1*c,funNum)/(i*c));
        g2 = exp(rq*a) * real(exp(-b0*i*d) * CZCharFun(S0,tau,a,params,K,rf,q,d,-b1*d,funNum)/(i*d));
        g3 = exp(rq*b) * real(exp(-b0*i*c) * CZCharFun(S0,tau,b,params,K,rf,q,c,-b1*c,funNum)/(i*c));
        g4 = exp(rq*b) * real(exp(-b0*i*d) * CZCharFun(S0,tau,b,params,K,rf,q,d,-b1*d,funNum)/(i*d));
        term1 = g1 + g2 + g3 + g4;
        h1 = exp(rq*(a+b)/2) * real(exp(-b0*i*c) * CZCharFun(S0,tau,(a+b)/2,params,K,rf,q,c,-b1*c,funNum)/(i*c));
        h2 = exp(rq*(a+b)/2) * real(exp(-b0*i*d) * CZCharFun(S0,tau,(a+b)/2,params,K,rf,q,d,-b1*d,funNum)/(i*d));
        h3 = exp(rq*a) * real(exp(-b0*i*(c+d)/2) * CZCharFun(S0,tau,a,params,K,rf,q,(c+d)/2,-b1*(c+d)/2,funNum)/(i*(c+d)/2));
        h4 = exp(rq*b) * real(exp(-b0*i*(c+d)/2) * CZCharFun(S0,tau,b,params,K,rf,q,(c+d)/2,-b1*(c+d)/2,funNum)/(i*(c+d)/2));
        term2 = h1 + h2 + h3 + h4;
        term3 = exp(rq*(a+b)/2) * real(exp(-b0*i*(c+d)/2) * CZCharFun(S0,tau,(a+b)/2,params,K,rf,q,(c+d)/2,-b1*(c+d)/2,funNum)/(i*(c+d)/2));
        Int(t,x) = (b-a)*(d-c)/16*(term1 + 2*term2 + 4*term3);
    end
end
y = sum(sum(Int));