function y = MNObjFun(param,param0,tau,tau0,S,K,rf,q,MktPrice,PutCall,MktIV,ObjFun,x,w)

% param    = Current parameter values to optimize
% param0   = Old parameter values
% tau      = Current maturity
% tau0     = Old maturity
% S        = Spot
% K        = Strike
% rf       = Risk Free Rate
% q        = Dividend yield
% MktPrice = Market Price of the option
% PutCall  = 'C'all or 'P'ut
% Weight   = Weights for the optimizer

% kappa  = param(1);
% theta  = param(2);
% sigma  = param(3);
% v0     = param(4);
% rho    = param(5);

% BlackScholes vega
BSV = @(S,K,r,q,v,T) (S*exp(-q*T)*normpdf((log(S/K) + T*(r+0.5*v^2))/v/sqrt(T))*sqrt(T));

[NK,NT] = size(MktPrice);
ModelPrice = zeros(NK,NT);
error      = zeros(NK,NT);

for k=1:NK
    MPrice(k) = MNPriceGaussLaguerre(param,param0,tau,tau0,K(k),S,PutCall(k,1),rf,q,x,w);
    switch ObjFun
        case 1
            % MSE
            error(k) =  (MPrice(k) - MktPrice(k))^2;
        case 2
            % RMSE
            error(k) =  (MPrice(k) - MktPrice(k))^2 / MktPrice(k);
        case 3
            % CHJ (2009)
            Vega(k) = BSV(S,K(k),rf,q,MktIV(k),tau);
            error(k) = sqrt((MPrice(k) - MktPrice(k))^2) / Vega(k);
    end
end

y = sum(sum(error)) / (NT*NK);
