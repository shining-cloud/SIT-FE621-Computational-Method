function y = HestonObjFun(param,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,ObjFun,Method,a,b,Tol,MaxIter,trap);

% S  = Spot price
% rf = risk free rate
% q  = dividend yield
% MktPrice = matrix of market prices.
% K = vector of strikes.
% T = vector of maturities.
% PutCall = matrix of put/call indicator 'C' or 'P'
% MktIV = matrix of market implied volatilities
% x = abscissas for Gauss Laguerre integration
% w = weights for Gauss Laguerre integration
% ObjFun = Type of Objective Function
%    1 = MSE
%    2 = Relative MSE
%    3 = Relative IVMSE
%    4 = Relative IVMSE Christoffersen, Heston, Jacobs proxy
% Method = Method used to obtain the option prices
%    1 = Heston with Gauss Laguerre integration
%    2 = Lewis vol-of-vol Series II expansion (IV directly)
%    3 = Lewis vol-of-vol Series I expansion
% a = Bisection algorithm, small initial estimate
% b = Bisection algorithm, large initial estimate
% Tol = Bisection algorithm, tolerance
% MaxIter = Bisection algorithm, maximum iterations
% trap = 1:"Little Trap" c.f., 0:original Heston c.f.

kappa  = param(1); 
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);
lambda = 0;

[NK,NT] = size(MktPrice);

for t=1:NT
	for k=1:NK
		% Select the method for obtaining the price
		if Method==1
			% Heston Call price with Gauss Laguerre integration
			CallPrice = HestonPriceGaussLaguerre('C',S,K(k),T(t),rf,q,param,trap,x,w);
		elseif Method==2
			% Lewis vol-of-vol Series II expansion Call price
			[iv CallPrice] = SeriesIICall(S,K(k),rf,q,T(t),v0,rho,theta,kappa,sigma,'C');
		elseif Method==3
			% Lewis vol-of-vol Series I expansion Call price
			[iv CallPrice] = SeriesICall(S,K(k),rf,q,T(t),v0,rho,theta,kappa,sigma);
		end
		% Obtain the call price or put price
		if PutCall(k,t)=='C'
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
			% IVMSE
			if Method==2 | Method==3
				ModelIV = iv;
			else
			ModelIV = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
			end
			error(k,t) = (ModelIV - MktIV(k,t))^2 / MktIV(k,t);
		elseif ObjFun == 4
			% IVRMSE Christoffersen, Heston, Jacobs proxy
			d = (log(S/K(k)) + (rf-q+MktIV(k,t)^2/2)*T(t))/MktIV(k,t)/sqrt(T(t));
			Vega(k,t) = S*normpdf(d)*sqrt(T(t));
			error(k,t) = (ModelPrice(k,t) - MktPrice(k,t))^2/Vega(k,t)^2;
		end
	end
end

y = sum(sum(error));


