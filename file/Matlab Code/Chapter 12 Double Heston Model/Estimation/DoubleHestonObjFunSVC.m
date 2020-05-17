function y = DoubleHestonObjFunSVC(param,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap)

% Objective function for Strike Vector Computation
% param = Heston parameters
% S  = Spot price
% rf = risk free rate
% q  = dividend yield
% MktPrice = vector of market prices.
% K = vector of strikes.
% T = vector of maturities.
% PutCall = vector of put/call indicator 'C' or 'P'
% MktIV = vector of market implied volatilities
% x = abscissas for Gauss Laguerre integration
% w = weights for Gauss Laguerre integration
% trap = 1 Little Trap 0  Heston formulation
% ObjFun = Type of Objective Function
%    1 = MSE
%    2 = RMSE
%    3 = IVMSE
%    4 = Christoffersen, Heston, Jacobs
% a = Bisection algorithm, small initial estimate
% b = Bisection algorithm, large initial estimate
% Tol = Bisection algorithm, tolerance
% MaxIter = Bisection algorithm, maximum iterations

[NK,NT] = size(MktPrice);

for t=1:NT
	for j=1:length(x)
		% Store the characteristic functions at each time step
		phi = x(j);
        f2(j) = DoubleHestonCF(phi  ,param,T(t),S,rf,q,trap);
        f1(j) = DoubleHestonCF(phi-i,param,T(t),S,rf,q,trap);
    end
	for k=1:NK
        for j=1:length(x);
            phi = x(j);
            int2(j) = w(j) * real(exp(-i*phi*log(K(k)))*f2(j)/i/phi);
            int1(j) = w(j) * real(exp(-i*phi*log(K(k)))*f1(j)/i/phi/S/exp((rf-q)*T(t)));
        end
        P1 = 1/2 + 1/pi*sum(int1);
        P2 = 1/2 + 1/pi*sum(int2);
        % The call price
        CallPrice = S*exp(-q*T(t))*P1 - K(k)*exp(-rf*T(t))*P2;
        if strcmp(PutCall(k,t),'C')
            ModelPrice(k,t) = CallPrice;
        else
			ModelPrice(k,t) = CallPrice - S*exp(-q*T(t)) + exp(-rf*T(t))*K(k);
		end
		% Select the objective function
		if ObjFun == 1
			% MSE
			error(k,t) = (MktPrice(k,t) - ModelPrice(k,t))^2;
		elseif ObjFun == 2
			% RMSE
			error(k,t) = (MktPrice(k,t) - ModelPrice(k,t))^2 / MktPrice(k,t);
		elseif ObjFun == 3
			% IVRMSE Christoffersen, Heston, Jacobs proxy
			d = (log(S/K(k)) + (rf-q+MktIV(k,t)^2/2)*T(t))/MktIV(k,t)/sqrt(T(t));
			BSVega(k,t) = S*normpdf(d)*sqrt(T(t));
			error(k,t) = (ModelPrice(k,t) - MktPrice(k,t))^2 / BSVega(k,t)^2;
		end
	end
clear f1 f2
end
y = sum(sum(error));

