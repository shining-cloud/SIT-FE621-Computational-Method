function y = HestonBGMObjFun(param,T,S,K,rf,q,MktPrice,MktIV,ObjFun,PutCall)

% param    = Current parameter values to optimize
% T        = Current maturity
% S        = Spot
% K        = Strike
% rf       = Risk Free Rate
% q        = Dividend yield
% MktPrice = Market Price of the option
% ObjFun   = Choice of loss function
% PutCall  = 'P'ut or 'C'all

% Parameter order
% [kappa v0 | theta(1) sigma(1) rho(1) | theta(2) sigma(2) rho(2)]

BSV = @(S,K,r,q,v,T) (S*exp(-q*T)*normpdf((log(S/K) + T*(r+0.5*v^2))/v/sqrt(T))*sqrt(T));

[NK,NT] = size(MktPrice);

for t=1:NT
	for k=1:NK
		MPrice(k,t) = BGMApproxPriceTD(param,T(1:t),S,K(k),rf,q,PutCall(k,t));
        switch ObjFun
            case 1
                error(k,t) =  (MPrice(k,t) - MktPrice(k,t))^2;
            case 2
                error(k,t) =  (MPrice(k,t) - MktPrice(k,t))^2 / MktPrice(k,t);
            case 3
                Vega(k,t) = BSV(S,K(k),rf,q,MktIV(k,t),T(t));
                error(k,t) = sqrt((MPrice(k,t) - MktPrice(k,t))^2) / Vega(k,t);
        end
	end
end

y = sum(sum(error));
