function y = BGMApproxPriceTD(param,T,S,K,rf,q,PutCall)

% param  = vector of parameters to estimate
% T  = vector of maturities
% S  = Spot price
% K  = Strike price
% rf = Risk free rate
% q  = Dividend yield
% PutCall  = 'P'ut or 'C'all

% Order of the parameters to estimate
% [kappa v0 | theta(1) sigma(1) rho(1) | theta(2) sigma(2) rho(2) | etc..]

% Create a vector of all the parameters
NT = length(T);
Nparam = 2 + 3*NT;
kappa = param(1);
v0    = param(2);
for i=3:3:Nparam
	j = floor(i/3);
	theta(j) = param(i);
	sigma(j) = param(i+1);
	rho(j)   = param(i+2);
end

% First set of coefficients
a1T(1) = -rho(1)*sigma(1)*(2*theta(1)*exp(kappa*T(1))+v0*T(1)*kappa+v0-theta(1)*T(1)*kappa*exp(kappa*T(1))-theta(1)*kappa*T(1)-2*theta(1)-v0*exp(kappa*T(1)))/exp(kappa*T(1))/kappa^2;
a2T(1) = -1/2*rho(1)^2*sigma(1)^2*(kappa^2*v0*T(1)^2+6*theta(1)*exp(kappa*T(1))+2*v0*T(1)*kappa+2*v0-kappa^2*theta(1)*T(1)^2-4*theta(1)*kappa*T(1)-2*theta(1)*T(1)*kappa*exp(kappa*T(1)) - 6*theta(1)-2*v0*exp(kappa*T(1)))/exp(kappa*T(1))/kappa^3;
b0T(1) = -1/4*sigma(1)^2*(5*theta(1)*exp(kappa*T(1))^2+4*exp(kappa*T(1))*v0*T(1)*kappa-4*theta(1)*T(1)*kappa*exp(kappa*T(1))-2*theta(1)*T(1)*kappa*exp(kappa*T(1))^2+2*v0 - theta(1)-2*v0*exp(kappa*T(1))^2-4*theta(1)*exp(kappa*T(1)))/exp(kappa*T(1))^2/kappa^3;
b2T(1) = a1T(1)^2/2;

A1 = 0; 
A2 = 0; 
B0 = 0;

% Remaining sets of coefficients
if NT>=2
    for t=1:NT-1
        % Coefficients for a1T(T+1)
        for j=1
            A1 = A1 + rho(j)*sigma(j)*(-exp(-kappa*T(t))+exp(-kappa*T(t+1)))*(theta(j)-v0*T(j)*kappa+theta(j)*T(j)*kappa-theta(j)*exp(kappa*T(j)))/kappa^2;
        end
        for j=2:t
            A1 = A1 + -rho(j)*sigma(j)*(-exp(-kappa*T(t))+exp(-kappa*T(t+1)))*(-v0*T(j-1)*kappa+theta(j)*T(j-1)*kappa-theta(j)*exp(kappa*T(j-1))+v0*T(j)*kappa-theta(j)*T(j)*kappa+theta(j)*exp(kappa*T(j)))/kappa^2;
        end
        for j=t+1
            A1 = A1 + -rho(j)*sigma(j)*(-v0*exp(kappa*T(t+1))-v0*T(t)*kappa*exp(kappa*T(t))+theta(j)*exp(kappa*T(t+1))+theta(j)*T(t)*kappa*exp(kappa*T(t))+theta(j)*T(t)*kappa*exp(kappa*(T(t)+T(t+1)))-theta(j)*exp(2*kappa*T(t))+exp(kappa*T(t))*v0-exp(kappa*(T(t)+T(t+1)))*theta(j)*T(t+1)*kappa+theta(j)*exp(kappa*(T(t)+T(t+1)))+exp(kappa*T(t))*v0*kappa*T(t+1)-theta(j)*exp(kappa*T(t))-exp(kappa*T(t))*theta(j)*kappa*T(t+1))/kappa^2*exp(-kappa*(T(t)+T(t+1)));
        end
        a1T(t+1) = a1T(t) + A1;
        
        % Coefficients for a2T(T+1)
        for j=1
            A2 = A2 -rho(j)^2*sigma(j)^2*(exp(-kappa*T(t))+exp(-kappa*T(t+1))*T(t)*kappa-exp(-kappa*T(t+1))-exp(-kappa*T(t+1))*kappa*T(t+1))*(theta(j)-v0*T(j)*kappa+theta(j)*T(j)*kappa-theta(j)*exp(kappa*T(j)))/kappa^3 ...
                  + 1/2*rho(j)^2*sigma(j)^2*(-exp(-kappa*T(t))+exp(-kappa*T(t+1)))*(2*kappa*theta(j)*T(t)+2*theta(j)-2*v0*T(t)*T(j)*kappa^2+v0*T(j)^2*kappa^2+2*theta(j)*T(t)*T(j)*kappa^2-theta(j)*T(j)^2*kappa^2-2*theta(j)*T(t)*exp(kappa*T(j))*kappa+2*theta(j)*exp(kappa*T(j))*kappa*T(j)-2*theta(j)*exp(kappa*T(j)))/kappa^3;
        end
        for j=2:t
            A2 = A2 + rho(j)^2*sigma(j)^2*(exp(-kappa*T(t))+exp(-kappa*T(t+1))*T(t)*kappa-exp(-kappa*T(t+1))-exp(-kappa*T(t+1))*kappa*T(t+1))*(-v0*T(j-1)*kappa+theta(j)*T(j-1)*kappa-theta(j)*exp(kappa*T(j-1))+v0*T(j)*kappa-theta(j)*T(j)*kappa+theta(j)*exp(kappa*T(j)))/kappa^3 ...
                  - 1/2*rho(j)^2*sigma(j)^2*(-exp(-kappa*T(t))+exp(-kappa*T(t+1)))*(-2*v0*T(t)*T(j-1)*kappa^2+v0*T(j-1)^2*kappa^2+2*theta(j)*T(t)*T(j-1)*kappa^2-theta(j)*T(j-1)^2*kappa^2-2*theta(j)*T(t)*exp(kappa*T(j-1))*kappa+2*theta(j)*exp(kappa*T(j-1))*kappa*T(j-1)-2*theta(j)*exp(kappa*T(j-1))+2*v0*T(t)*T(j)*kappa^2-v0*T(j)^2*kappa^2-2*theta(j)*T(t)*T(j)*kappa^2+theta(j)*T(j)^2*kappa^2+2*theta(j)*T(t)*exp(kappa*T(j))*kappa-2*theta(j)*exp(kappa*T(j))*kappa*T(j)+2*theta(j)*exp(kappa*T(j)))/kappa^3;
        end
        for j=t+1
            A2 = A2 + 1/2*rho(j)^2*sigma(j)^2*(2*v0*exp(kappa*T(t+1))-v0*T(t)^2*kappa^2*exp(kappa*T(t))+2*v0*T(t)*kappa*exp(kappa*T(t))+2*v0*kappa^2*T(t+1)*T(t)*exp(kappa*T(t))-2*theta(j)*exp(kappa*T(t+1))+theta(j)*T(t)^2*kappa^2*exp(kappa*T(t))-2*theta(j)*T(t)*kappa*exp(kappa*T(t))-2*theta(j)*kappa^2*T(t+1)*T(t)*exp(kappa*T(t))-2*theta(j)*T(t)*kappa*exp(kappa*(T(t)+T(t+1)))-2*theta(j)*exp(2*kappa*T(t))*kappa*T(t)+4*theta(j)*exp(2*kappa*T(t))+2*theta(j)*T(t+1)*exp(2*kappa*T(t))*kappa-2*exp(kappa*T(t))*v0-2*exp(kappa*T(t))*v0*kappa*T(t+1)+2*exp(kappa*(T(t)+T(t+1)))*theta(j)*T(t+1)*kappa-exp(kappa*T(t))*kappa^2*v0*T(t+1)^2-4*theta(j)*exp(kappa*(T(t)+T(t+1)))+exp(kappa*T(t))*kappa^2*theta(j)*T(t+1)^2+2*theta(j)*exp(kappa*T(t))+2*exp(kappa*T(t))*theta(j)*kappa*T(t+1))/kappa^3*exp(-kappa*(T(t)+T(t+1)));
        end
        a2T(t+1) = a2T(t) + A2;
        
        % Coefficients for b0T(T+1)
        for j=1
            B0 = B0 -1/4*sigma(j)^2*(-exp(-2*kappa*T(t))+2*exp(-kappa*(T(t)+T(t+1)))-exp(-2*kappa*T(t+1)))*(-2*v0+theta(j)+2*v0*exp(kappa*T(j))-2*theta(j)*exp(kappa*T(j))+theta(j)*exp(2*kappa*T(j)))/kappa^3 ...
                  + 1/2*sigma(j)^2*(-theta(j)*exp(kappa*T(t+1))+theta(j)*exp(kappa*T(t))+2*v0*exp(kappa*T(t+1))-2*exp(kappa*T(t))*v0-2*theta(j)*exp(kappa*(T(t)+T(t+1)))+2*theta(j)*exp(2*kappa*T(t))+2*exp(kappa*(T(t)+T(t+1)))*v0*T(j)*kappa-2*exp(kappa*(T(t+1)+T(j)))*v0-2*exp(kappa*(T(t)+T(t+1)))*theta(j)*T(j)*kappa+2*exp(kappa*(T(t+1)+T(j)))*theta(j)+2*exp(kappa*(T(t)+T(t+1)+T(j)))*theta(j)-exp(kappa*(T(t+1)+2*T(j)))*theta(j)-2*exp(2*kappa*T(t))*v0*T(j)*kappa+2*exp(kappa*(T(t)+T(j)))*v0+2*exp(2*kappa*T(t))*theta(j)*T(j)*kappa-2*exp(kappa*(T(t)+T(j)))*theta(j)-2*exp(kappa*(2*T(t)+T(j)))*theta(j)+exp(kappa*(T(t)+2*T(j)))*theta(j))*exp(-kappa*(2*T(t)+T(t+1)))/kappa^3;
        end
        for j=2:t
            B0 = B0 + -1/4*sigma(j)^2*(-exp(-2*kappa*T(t))+2*exp(-kappa*(T(t)+T(t+1)))-exp(-2*kappa*T(t+1)))*(-2*v0*exp(kappa*T(j-1))+2*theta(j)*exp(kappa*T(j-1))-theta(j)*exp(2*kappa*T(j-1))+2*v0*exp(kappa*T(j))-2*theta(j)*exp(kappa*T(j))+theta(j)*exp(2*kappa*T(j)))/kappa^3 ...
                  + 1/2*sigma(j)^2*(exp(-kappa*T(t))-exp(-kappa*T(t+1)))*(-2*v0*T(j-1)*kappa*exp(kappa*T(t))+2*v0*exp(kappa*T(j-1))+2*theta(j)*T(j-1)*kappa*exp(kappa*T(t))-2*theta(j)*exp(kappa*T(j-1))-2*theta(j)*exp(kappa*(T(j-1)+T(t)))+theta(j)*exp(2*kappa*T(j-1))+2*v0*T(j)*kappa*exp(kappa*T(t))-2*v0*exp(kappa*T(j))-2*theta(j)*T(j)*kappa*exp(kappa*T(t))+2*theta(j)*exp(kappa*T(j))+2*exp(kappa*(T(t)+T(j)))*theta(j)-theta(j)*exp(2*kappa*T(j)))*exp(-kappa*T(t))/kappa^3;
        end
        for j=t+1
            B0 = B0 + -1/4*sigma(j)^2*(-2*v0*exp(2*kappa*T(t+1))-4*v0*T(t)*kappa*exp(kappa*(T(t)+T(t+1)))+2*v0*exp(2*kappa*T(t))+2*theta(j)*exp(2*kappa*T(t+1))+4*theta(j)*T(t)*kappa*exp(kappa*(T(t)+T(t+1)))-2*theta(j)*exp(2*kappa*T(t))+2*theta(j)*T(t)*kappa*exp(kappa*(T(t)+2*T(t+1)))-4*theta(j)*exp(kappa*(2*T(t)+T(t+1)))+theta(j)*exp(3*kappa*T(t))+3*theta(j)*exp(kappa*(T(t)+2*T(t+1)))+4*exp(kappa*(T(t)+T(t+1)))*v0*kappa*T(t+1)-4*exp(kappa*(T(t)+T(t+1)))*theta(j)*T(t+1)*kappa-2*exp(kappa*(T(t)+2*T(t+1)))*theta(j)*T(t+1)*kappa)/kappa^3*exp(-kappa*(T(t)+2*T(t+1)));
        end
        b0T(t+1) = b0T(t) +(B0);

        % Coefficients for b2T(t+1)
        b2T(t+1) = a1T(t+1)^2/2;
    end
end

% Coefficients for the expansion are the last ones in the iterations
A1T = a1T(end);
A2T = a2T(end);
B0T = b0T(end);
B2T = b2T(end);

% log Spot price
x = log(S);

% Integrated variance
wT = (v0-theta(end))*(1-exp(-kappa*T(end)))/kappa + theta(end)*T(end);
y = wT;

% Black Scholes Put Price
g = y^(-1/2) * (-x + log(K) - (rf-q)*T(end)) - (1/2)*sqrt(y);
f = y^(-1/2) * (-x + log(K) - (rf-q)*T(end)) + (1/2)*sqrt(y);
BSPut = K*exp(-rf*T(end))*normcdf(f) - S*exp(-q*T(end))*normcdf(g);

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
PHIgx   = -phig*y^(-1/2);
PHIgy   = -(1/2)*f*phig/y;
PHIgxy  =  (1/2)*y^(-3/2)*phig*(1-f*g);
PHIgx2y =  (1/2)*y^(-2)*phig*(2*g+f-g^2*f);
PHIgy2  =  (1/2)*y^(-2)*phig*(f+g/2-g*f^2/2);
PHIgxy2 =  (1/2)*y^(-2)*(phigx*(f+g/2-f^2*g/2) + phig*(fx+gx/2-f*fx*g/2-f^2*gx/2));
PHIgx2y2 = (1/2)*((y^(-2)*phigy-2*y^(-3)*phig)*(2*g+f-g^2*f) +...
                   y^(-2)*phig*(2*gy+fy-2*g*gy*f-g^2*fy));

% Derivatives of Black-Scholes Put
dPdxdy   = K*exp(-rf*T(end))*PHIfxy   - exp(-q*T(end))*S*(PHIgy + PHIgxy);
dPdx2dy  = K*exp(-rf*T(end))*PHIfx2y  - exp(-q*T(end))*S*(PHIgy + 2*PHIgxy + PHIgx2y);
dPdy2    = K*exp(-rf*T(end))*PHIfy2   - exp(-q*T(end))*S*PHIgy2;
dPdx2dy2 = K*exp(-rf*T(end))*PHIfx2y2 - exp(-q*T(end))*S*(PHIgy + 2*PHIgxy + PHIgx2y + PHIgy2 + 2*PHIgxy2 + PHIgx2y2);

% Benhamou, Gobet, Miri expansion
Put = BSPut + A1T*dPdxdy + A2T*dPdx2dy + B0T*dPdy2 + B2T*dPdx2dy2;

% Return the put or the call by put-call parity
if strcmp(PutCall(1),'P')
	y = Put;
elseif strcmp(PutCall(1),'C')
	y = Put - K*exp(-rf*T(end)) + S*exp(-q*T(end));
end

