function [Price y theta] = MSPriceBS(S,K,T,sigma,r,q)

% Moneyness
theta = log(K/S)/sigma/sqrt(T);

% Optimization settings
start = 2;
options = optimset('LargeScale','off');

% Find the barrier
[y feval] = fmincon(@(p) -MSPutBS(p,theta,K,sigma,r,q,T),start,[],[],[],[],theta,[],[],options);

% American put Black Scholes price
Price = MSPutBS(y,theta,K,sigma,r,q,T);
