function y = HestonObjFunFRFT(param,S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,...
	         ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule)

% param = Heston parameter vector         
% S  = Spot price
% K1 = starting value for strike grid
% rf = risk free rate
% q  = dividend yield
% MktPrice = matrix of market prices.
% K = vector of strikes.
% T = vector of maturities.
% PutCall = matrix of put/call indicator 'C' or 'P'
% MktIV = matrix of market implied volatilities
% ObjFun = Type of Objective Function
%    1 = MSE
%    2 = Relative MSE
%    3 = Relative IVMSE
%    4 = Relative IVMSE Christoffersen, Heston, Jacobs proxy
% a = Bisection algorithm, small initial estimate
% b = Bisection algorithm, large initial estimate
% Tol = Bisection algorithm, tolerance
% MaxIter = Bisection algorithm, maximum iterations
% trap = 1:"Little Trap" c.f., 0:original Heston c.f.
% FFT settings
%   N  = number of discretization points
%   uplimit = Upper limit of integration
%   eta = integration grid size
%   rule = integration rule
%     rule = 'T' --> Trapezoidal
%     rule = 'S' --> Simpson's Rule

kappa  = param(1); 
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);
lambda = 0;

[NK,NT] = size(MktPrice);
lambdainc = 2/N*log(S/K1);

for t=1:NT
    [CallFRFT KK lambdainc eta] = HestonCallFRFT(N,uplimit,S,rf,q,T(t),kappa,theta,lambda,rho,sigma,v0,trap,alpha,rule,eta,lambdainc);
%	CallPrice = interp1(KK,CallFRFT,K,'linear')';
    CallPrice = LinearInterpolate(KK,CallFRFT,K);
    for k=1:NK
		if strcmp(PutCall(k,t),'C')
			ModelPrice(k,t) = CallPrice(k);
		else
			ModelPrice(k,t) = CallPrice(k) - S*exp(-q*T(t)) + exp(-rf*T(t))*K(k);
		end
	end
	% Select the objective function
	if ObjFun == 1
		% MSE
		error(:,t) = (MktPrice(:,t) - ModelPrice(:,t)).^2;
	elseif ObjFun == 2
		% RMSE
		error(:,t) = (MktPrice(:,t) - ModelPrice(:,t)).^2 ./ MktPrice(:,t);
	elseif ObjFun == 3
		for k=1:NK
			% IVMSE
            ModelIV = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
			error(k,t) = (ModelIV - MktIV(k,t))^2 / MktIV(k,t);
		end
	elseif ObjFun == 4
		for k=1:NK
			% IVRMSE Christoffersen, Heston, Jacobs proxy
			d = (log(S/K(k)) + (rf-q+MktIV(k,t)^2/2)*T(t))/MktIV(k,t)/sqrt(T(t));
			Vega(k,t) = S*normpdf(d)*sqrt(T(t));
			error(k,t) = (ModelPrice(k,t) - MktPrice(k,t))^2/Vega(k,t)^2;
		end
 	end
end

y = sum(sum(error)) / (NT*NK);


