function [S V SimPrice] = QEPrice(params,PutCall,S0,K,Mat,r,q,T,N,gamma1,gamma2,MC,phic,icdf)

% Heston Call or Put price using Kahl-Jackel IJK or Pathwise schemes
% INPUTS
%   params = Heston parameters
%   PutCall = 'C' for Call option, 'P' for put
%   S0  = Spot Price
%   K = Strike price
%   Mat = Maturity
%   r = risk free rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
%   gamma1 = first gamma parameter
%   gamma2 = second gamma parameter
%   MC = 1:Martingale correction, 0:No correction
%   phic = cut-off for low and high regimes for v(t)
%   icdf = 1=Matlab algorithm
%          2=Wichura algorithm
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances
%   SimPrice = option price by simulation

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);

% Obtain the simulated stock price and simulated variance
[S V] = QESim(params,gamma1,gamma2,S0,Mat,r,q,T,N,MC,phic,icdf);

% Terminal stock prices
ST = S(end,:);

% Payoff vectors
if strcmp(PutCall,'C')
	Payoff = max(ST - K,0);
elseif strcmp(PutCall,'P')
	Payoff = max(K - ST,0);
end
	
% Simulated price
SimPrice = exp(-r*Mat)*mean(Payoff);
