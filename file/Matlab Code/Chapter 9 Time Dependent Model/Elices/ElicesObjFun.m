function y = ElicesObjFun(param,paramfixed,v0,T,S,K,rf,q,MktPrice,MktIV,PutCall,ObjFun,x,w,trap)

% param    = Current period parameter values to optimize
% paramfixed = Previous period parameter values
% T        = Maturities (shortest first) T(1),T(2),...,T(N-1),T(N)
% S        = Spot
% K        = Strike
% rf       = Risk Free Rate
% q        = Dividend yield
% MktPrice = Market Price of the option
% MktIV    = Market implied vol of the option
% PutCall  = 'C'all or 'P'ut
% ObjFun   = Choice of objective function
% x = integration abscissas
% w = integration weights

% Maturities are in a vector, shortest maturities first
%    T = (T(1),T(2),...., T(N-1),T(N))
% Fixed parameters are in a matrix, those for the shortest maturity first
%    paramfixed: (paramfixed(1),paramfixed(2), ...,paramfixed(N-1))
% Dimension of parameter rows is one less than dimension of maturities

[NK,NT] = size(MktPrice);

% BlackScholes vega
BSV = @(S,K,r,q,v,T) (S*exp(-q*T)*normpdf((log(S/K) + T*(r+0.5*v^2))/v/sqrt(T))*sqrt(T));

for k=1:NK
    MPrice(k) = ElicesPrice(PutCall,S,K(k),T,rf,q,param,paramfixed,v0,trap,x,w);
    switch ObjFun
        case 1
            % MSE
            error(k) = (MktPrice(k) - MPrice(k))^2;
        case 2
            % RMSE
            error(k) = (MktPrice(k) - MPrice(k))^2 / MktPrice(k);
        case 3
            % CHJ (2009)
            Vega(k) = BSV(S,K(k),rf,q,MktIV(k),T(end));
            error(k) = sqrt((MPrice(k) - MktPrice(k))^2) / Vega(k);
    end
end

y = sum(sum(error));



