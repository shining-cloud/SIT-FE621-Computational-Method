function AmerPut = MSPrice(S,Strike,T,r,q,params,trap,method,A,B,N,NumTerms,yinf)

% Function to calculate American put price by control variate method
% with Medvedev-Scaillet expansion
% MSPut = function for calculating MS put price by expansion
% Option features
%    S = spot price
%    Strike = strike price
%    T = maturity
%    r = risk free rate
%    q = dividend yield
% Heston parameters and settings
%    kappaV = mean reversion rate
%    thetaV = mean reversion level
%    lambda = volatility risk
%    v0     = initial variance
%    rho    = correlation
%    trap   = "Little Trap" formulation
% Newton-Cotes integration for European put
%    method = Integration method
%    A = Lower limit
%    B = Upper limit
%    N = Number of points
% Settings for MS expansion
%    NumTerms = number of terms (3,4, or 5)
%    yinf = y infinity for European put

% Parameter vector
kappaV = params(1);
thetaV = params(2);
sigmaV = params(3);
v0     = params(4);
rho    = params(5);
lambda = 0;

% Closed-form European put
EuroPutClosed = HestonPriceNewtonCoates('P',S,Strike,T,r,q,kappaV,thetaV,sigmaV,lambda,v0,rho,trap,method,A,B,N);

% Moneyness
theta = log(Strike/S)/sqrt(v0)/sqrt(T);

% Find the barrier level, y
lo = max(2,theta);
y = fminbnd(@(p) -MSPutHeston(p,theta,Strike,params,r,q,T,NumTerms),lo,4);
if y<theta
    y = theta;
end

% Euro and Amer put by MS expansion
EuroPutMS = MSPutHeston(yinf,theta,Strike,params,r,q,T,NumTerms);
AmerPutMS = MSPutHeston(y,   theta,Strike,params,r,q,T,NumTerms);

% Find the early exercise premium
EEP = AmerPutMS - EuroPutMS;

% American put terms
AmerPut = EuroPutClosed + EEP;
