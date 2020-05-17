function y = BGMApproxPrice(params,S,K,rf,q,T,trap,PutCall)

kappa  = params(1);
theta  = params(2);
sigma  = params(3);
rho    = params(4);
v0     = params(5);
lambda = params(6);

% log Spot price
x = log(S);

% Integrated variance
wT = (v0-theta)*(1-exp(-kappa*T))/kappa + theta*T;
y = wT;

% Black Scholes Put Price
g = y^(-1/2) * (-x + log(K) - (rf-q)*T) - (1/2)*sqrt(y);
f = y^(-1/2) * (-x + log(K) - (rf-q)*T) + (1/2)*sqrt(y);
BSPut = K*exp(-rf*T)*normcdf(f) - S*exp(-q*T)*normcdf(g);

% Shortcut notation
k  = kappa;
kT = kappa*T;
ekT  = exp( k*T);
ekTm = exp(-k*T);

% Coefficients for the expansion
a1T = (rho*sigma*ekTm/k^2) * (v0*(-kT+ekT-1) + theta*(kT+ekT*(kT-2)+2));
a2T = (rho^2*sigma^2*ekTm/2/k^3) * (v0*(-kT*(kT+2)+2*ekT-2) + theta*(2*ekT*(kT-3)+kT*(kT+4)+6));
b0T = (sigma^2*exp(-2*kT)/4/k^3) * (v0*(-4*ekT*kT+2*exp(2*kT)-2) + theta*(4*ekT*(kT+1)+exp(2*kT)*(2*kT-5)+1));
b2T = a1T^2/2;

% Normal pdf, phi(f) and phi(g)
phif = exp(-f^2/2)/sqrt(2*pi);
phig = exp(-g^2/2)/sqrt(2*pi);

% Derivatives of f and g
fx = -y^(-1/2);
fy = -1/2/y*g;
gx = fx;
gy = -1/2/y*f;

% The cdf PHI(f) and PHI(g)
PHIf = normcdf(f);
PHIg = normcdf(g);

% Derivatives of the pdf phi(f)
phifx = y^(-1/2)*f*phif;
phify = (1/2)/y*f*g*phif;

% Derivatives of the cdf PHI(f)
PHIfxy   = (1/2)*y^(-3/2)*phif*(1-f*g);
PHIfx2y  = (1/2)*y^(-2)*phif*(2*f+g-f^2*g);
PHIfy2   = (1/2)*y^(-2)*phif*(g+f/2-f*g^2/2);
PHIfx2y2 = (1/2)*((y^(-2)*phify-2*y^(-3)*phif)*(2*f+g-f^2*g) +...
	               y^(-2)*phif*(2*fy+gy-2*f*fy*g-f^2*gy));
			   
% Derivatives of the pdf phi(g)
phigx = y^(-1/2)*g*phig;
phigy = (1/2)/y*f*g*phig;

% Derivatives of cdf PHI(g)
PHIgx = -phig*y^(-1/2);
PHIgy   = -(1/2)*f*phig/y;
PHIgxy  =  (1/2)*y^(-3/2)*phig*(1-f*g);
PHIgx2y =  (1/2)*y^(-2)*phig*(2*g+f-g^2*f);
PHIgy2  =  (1/2)*y^(-2)*phig*(f+g/2-g*f^2/2);
PHIgxy2 =  (1/2)*y^(-2)*(phigx*(f+g/2-f^2*g/2) + phig*(fx+gx/2-f*fx*g/2-f^2*gx/2));
PHIgx2y2 = (1/2)*((y^(-2)*phigy-2*y^(-3)*phig)*(2*g+f-g^2*f) +...
                   y^(-2)*phig*(2*gy+fy-2*g*gy*f-g^2*fy));

% Derivatives of Black-Scholes Put
dPdxdy   = K*exp(-rf*T)*PHIfxy   - exp(-q*T)*S*(PHIgy + PHIgxy);
dPdx2dy  = K*exp(-rf*T)*PHIfx2y  - exp(-q*T)*S*(PHIgy + 2*PHIgxy + PHIgx2y);
dPdy2    = K*exp(-rf*T)*PHIfy2   - exp(-q*T)*S*PHIgy2;
dPdx2dy2 = K*exp(-rf*T)*PHIfx2y2 - exp(-q*T)*S*(PHIgy + 2*PHIgxy + PHIgx2y + PHIgy2 + 2*PHIgxy2 + PHIgx2y2);

% Benhamou, Gobet, Miri expansion
Put = BSPut + a1T*dPdxdy + a2T*dPdx2dy + b0T*dPdy2 + b2T*dPdx2dy2;

% Return the put or the call by put-call parity
if strcmp(PutCall(1),'P')
	y = Put;
else
	y = Put - K*exp(-rf*T) + S*exp(-q*T);
end


