function y = DoubleGaussLegendre(S0,tau,params,K,rf,q,b0,b1,xt,wt,xs,ws,a,b,c,d,funNum)

% Double integration by composite Gauss Legendre rule
%   S0  = spot price
%   tau = maturity
%   params = Heston parameters
%   K   = strike price
%   rf  = risk free rate
%   q   = dividend yield
%   b0 = Chiarella and Ziogas parameter
%   b1 = Chiarella and Ziogas parameter
%   xs = abscissas for Gauss-Laguerre
%   ws = weights for Gauss-Laguerre
%   xt = abscissas for Gauss-Legendre
%   wt = weights for Gauss-Legendre
%   integration grid (t,x) on (a,b) x (c,d)
%   funNum = c.f. 1 or c.f. 2

Nt = length(xt);
Ns = length(xs);
h1 = (b-a)/2;
h2 = (b+a)/2;
k1 = (d-c)/2;
k2 = (d+c)/2;
y = 0;

if (funNum == 1)
    qr = q;
elseif (funNum == 2)
    qr = rf;
end

for t=1:Nt;
    time = h1*xt(t) + h2;
    for x=1:Ns
        phi = k1*xs(x) + k2;
        fun = exp(qr*time) * real(exp(-b0*i*phi) * CZCharFun(S0,tau,time,params,K,rf,q,phi,-b1*phi,funNum)/(i*phi));
        y = y + h1*k1*wt(t)*ws(x) * fun;
    end
end


