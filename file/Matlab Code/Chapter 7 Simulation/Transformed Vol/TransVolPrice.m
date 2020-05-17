function [S v SimPrice] = TransVolPrice(scheme,params,PutCall,S0,K,Mat,r,q,T,N)

% Heston Call or Put price using
% Zhu (2010) Transformed Volatility discretization of the Heston model
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
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances
%   F = Matrix of overriden negative variances
%   SimPrice = option price by simulation

% Obtain the simulated stock price and simulated variance
[S v] = TransVolSim(scheme,params,S0,Mat,r,q,T,N);

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
