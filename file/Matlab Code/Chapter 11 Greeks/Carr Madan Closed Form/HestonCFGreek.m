function df = HestonCFGreek(phi,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Trap,Greek)

% Returns the Heston characteristic function, f2 (the second CF)
% Uses original Heston formulation 1, or formulation 2 from "The Little Heston Trap"
% phi = integration variable
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0     = initial variance
% Option features.
%    PutCall = 'C'all or 'P'ut
%    K = strike price
%    S = spot price
%    r = risk free rate
% Trap = 0 --> Original Heston Formulation
% Trap = 1 --> Little Trap formulation

% Log of the stock price.
x0 = log(S);

% Parameter a
a = kappa.*theta;

% Parameters u, b, d, and b
u = -0.5;
b = kappa + lambda;
d = sqrt((rho.*sigma.*i.*phi - b).^2 - sigma.^2.*(2.*u.*i.*phi - phi.^2));
g = (b - rho.*sigma.*i.*phi + d) ./ (b - rho.*sigma.*i.*phi - d);

if Trap==1
	% "Little Heston Trap" formulation
	c = 1./g;
	G = (1 - c.*exp(-d.*tau))./(1-c);
	C = r.*i.*phi.*tau + a./sigma.^2.*((b - rho.*sigma.*i.*phi - d).*tau - 2.*log(G));
	D = (b - rho.*sigma.*i.*phi - d)./sigma.^2.*((1-exp(-d.*tau))./(1-c.*exp(-d.*tau)));
elseif Trap==0
	% Original Heston formulation.
	G = (1 - g.*exp(d.*tau))./(1-g);
	C = r.*i.*phi.*tau + a./sigma.^2.*((b - rho.*sigma.*i.*phi + d).*tau - 2.*log(G));
	D = (b - rho.*sigma.*i.*phi + d)./sigma.^2.*((1-exp(d.*tau))./(1-g.*exp(d.*tau)));
end

% The characteristic function.
f = exp(C + D.*v0 + i.*phi.*x0);

% Derivatives of f2 required for the Greeks
if strcmp(Greek,'Price');
    df = f;
elseif strcmp(Greek,'Delta')
    df = f*i*phi/S;
elseif strcmp(Greek,'Gamma')
    df = i*phi*f*(i*phi - 1)/S^2;
elseif strcmp(Greek,'Theta')
	dDdT = d*exp(d*tau)*(b-rho*sigma*phi*i+d)*(g-1)/sigma^2/(1-g*exp(d*tau))^2;
	dCdT = r*phi*i + kappa*theta/sigma^2 ...
	     * ((b-rho*sigma*phi*i+d) + 2*g*d*exp(d*tau)/(1-g*exp(d*tau)));
	df = r*f - f*(dCdT + dDdT*v0);
elseif strcmp(Greek,'Rho')
    dCdr = i*phi*tau;
    df = (f*dCdr - tau*f);
elseif strcmp(Greek,'Vega1')
    df = f*D;
elseif strcmp(Greek,'Vanna')
    df = f*i*phi*D/S;
elseif strcmp(Greek,'Volga')
    df = 2*D*f*(2*D*v0 + 1);
end

