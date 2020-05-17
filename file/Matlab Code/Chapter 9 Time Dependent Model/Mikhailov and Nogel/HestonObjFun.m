function y = HestonObjFun(param,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap)

% BlackScholes vega
BSV = @(S,K,r,q,v,T) (S*exp(-q*T)*normpdf((log(S/K) + T*(r+0.5*v^2))/v/sqrt(T))*sqrt(T));

[NK,NT] = size(MktPrice);
MPrice = zeros(NK,NT);
error  = zeros(NK,NT);
Vega   = zeros(NK,NT);
for t=1:NT
	for k=1:NK
        MPrice(k,t) = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,param,trap,x,w);
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

y = sum(sum(error)) / (NT*NK);




