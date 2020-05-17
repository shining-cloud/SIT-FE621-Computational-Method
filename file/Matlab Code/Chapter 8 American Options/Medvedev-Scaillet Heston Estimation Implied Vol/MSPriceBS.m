function Price = MSPriceBS(S,K,T,sigma,r,q)

% Log moneyness
theta = log(K/S)/sigma/sqrt(T);

% Starting value
start = 2;
options = optimset('Display','off','Algorithm','active-set');

% Find the barrier using constrained optimization
[y feval] = fmincon(@(p) -MSPutBS(p,theta,K,sigma,r,q,T),start,[],[],[],[],theta,[],[],options);

% The American put price under Black Scholes
Price = MSPutBS(y,theta,K,sigma,r,q,T);
