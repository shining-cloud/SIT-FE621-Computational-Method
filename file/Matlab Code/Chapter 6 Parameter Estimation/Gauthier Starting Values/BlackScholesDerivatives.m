function [dPdxdy dPdx2dy dPdy2 dPdx2dy2] = BlackScholesDerivatives(kappa,theta,v0,S,K,T,rf,q)

% dPdxdy   = (1,1)
% dPdx2dy  = (2,1)
% dPdy2    = (0,2)
% dPdx2dy2 = (2,2)

% log Spot price
x = log(S);

% Integrated variance
wT = (v0-theta)*(1-exp(-kappa*T))/kappa + theta*T;
y = wT;

% Black Scholes Put Price
g = y^(-1/2)*(-x + log(K) - (rf-q)*T) - (1/2)*sqrt(y);
f = y^(-1/2)*(-x + log(K) - (rf-q)*T) + (1/2)*sqrt(y);
BSPut = K*exp(-rf*T)*normcdf(f) - S*exp(-q*T)*normcdf(g);

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

% Derivatives obtained using the symbolic toolbox in Matlab
% dPdxdy   = -1/2*K*exp(-rf*T)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)+1/4*K*exp(-rf*T)/pi^(1/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(3/2)+1/2*S*exp(-T*q)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)-1/4*S*exp(-T*q)/pi^(1/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(3/2);
% dPdx2dy  = -1/2*K*exp(-rf*T)/pi^(1/2)*(1/2/y-1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)-1/4*K*exp(-rf*T)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+1/8*K*exp(-rf*T)/pi^(1/2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(5/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+1/2*S*exp(-T*q)/pi^(1/2)*(1/2/y-1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)+1/4*S*exp(-T*q)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)-1/8*S*exp(-T*q)/pi^(1/2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(5/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2);
% dPdy2    =  K*exp(-rf*T)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*(1/4*2^(1/2)/y^(1/2)-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)*2^(1/2)/y^(3/2))+K*exp(-rf*T)/pi^(1/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*(-1/4*2^(1/2)/y^(3/2)+3/16*(-2*x+2*log(K)-2*rf*T+2*T*q+y)*2^(1/2)/y^(5/2))-S*exp(-T*q)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*(1/4*2^(1/2)/y^(1/2)-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)*2^(1/2)/y^(3/2))-S*exp(-T*q)/pi^(1/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*(-1/4*2^(1/2)/y^(3/2)+3/16*(-2*x+2*log(K)-2*rf*T+2*T*q+y)*2^(1/2)/y^(5/2));
% dPdx2dy2 = -1/2*K*exp(-rf*T)/pi^(1/2)*(-1/y^2+(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^3)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)-1/2*K*exp(-rf*T)/pi^(1/2)*(1/2/y-1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)+1/4*K*exp(-rf*T)/pi^(1/2)*(1/2/y-1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(3/2)-1/4*K*exp(-rf*T)/pi^(1/2)*(-1/4/y+1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^3)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)-1/4*K*exp(-rf*T)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+1/2*K*exp(-rf*T)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(5/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)-1/4*K*exp(-rf*T)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)^2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+1/8*K*exp(-rf*T)/pi^(1/2)/y^(5/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)-5/16*K*exp(-rf*T)/pi^(1/2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(7/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+1/2*S*exp(-T*q)/pi^(1/2)*(-1/y^2+(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^3)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)+1/2*S*exp(-T*q)/pi^(1/2)*(1/2/y-1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(1/2)-1/4*S*exp(-T*q)/pi^(1/2)*(1/2/y-1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)/y^(3/2)+1/4*S*exp(-T*q)/pi^(1/2)*(-1/4/y+1/2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^2-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^3)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+1/4*S*exp(-T*q)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)-1/2*S*exp(-T*q)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(5/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+1/4*S*exp(-T*q)/pi^(1/2)*(-1/4*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y+1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y^2)^2*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(3/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)-1/8*S*exp(-T*q)/pi^(1/2)/y^(5/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2)+5/16*S*exp(-T*q)/pi^(1/2)*(-2*x+2*log(K)-2*rf*T+2*T*q+y)/y^(7/2)*exp(-1/8*(-2*x+2*log(K)-2*rf*T+2*T*q+y)^2/y)*2^(1/2);
  

