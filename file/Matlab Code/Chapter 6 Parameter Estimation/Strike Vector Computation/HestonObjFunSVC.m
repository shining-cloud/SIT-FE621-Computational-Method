function y = HestonObjFunSVC(param,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,CF)

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
% CF = 'Heston' or 'Attari' characteristic function

kappa  = param(1); 
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);
lambda = 0;

[NK,NT] = size(MktPrice);

for t=1:NT
	for j=1:length(x)
		% Store the characteristic functions at each time step
		phi = x(j);
        if strcmp(CF,'Heston')
            f2(j) = HestonCF(phi,  kappa,theta,lambda,rho,sigma,T(t),S,rf,q,v0,trap);
            f1(j) = HestonCF(phi-i,kappa,theta,lambda,rho,sigma,T(t),S,rf,q,v0,trap) ...
                / (S*exp((rf-q)*T(t)));
        elseif strcmp(CF,'Attari')
            f(j) = AttariCF(phi,kappa,theta,lambda,rho,sigma,T(t),S,rf,q,v0,trap);
        end
    end
	for k=1:NK
        L = log(exp(-rf*T(t))*K(k)/S);
 		for j=1:length(x);
 			phi = x(j);        
            if strcmp(CF,'Heston')
                int1(j) = w(j) * real(exp(-i*phi*log(K(k)))*f1(j)/i/phi);
                int2(j) = w(j) * real(exp(-i*phi*log(K(k)))*f2(j)/i/phi);
            elseif strcmp(CF,'Attari')
                int1(j) = w(j) * ((real(f(j)) + imag(f(j))/phi)*cos(L*phi) ...
                        + (imag(f(j)) - real(f(j))/phi)*sin(L*phi)) / (1+phi^2);
            end
        end
        % The call price
        if strcmp(CF,'Heston')
            P1 = 1/2 + 1/pi*sum(int1);
            P2 = 1/2 + 1/pi*sum(int2);
            CallPrice = S*exp(-q*T(t))*P1 - K(k)*exp(-rf*T(t))*P2;
        elseif strcmp(CF,'Attari')
            CallPrice = S*exp(-q*T(t)) - K(k)*exp(-rf*T(t))*(1/2 + 1/pi*sum(int1));
        end
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
			% IVMSE
			ModelIV = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
			error(k,t) = (ModelIV - MktIV(k,t))^2;
		elseif ObjFun == 4
			% IVRMSE Christoffersen, Heston, Jacobs proxy
			d = (log(S/K(k)) + (rf-q+MktIV(k,t)^2/2)*T(t))/MktIV(k,t)/sqrt(T(t));
			Vega(k,t) = S*normpdf(d)*sqrt(T(t));
			error(k,t) = (ModelPrice(k,t) - MktPrice(k,t))^2 / Vega(k,t)^2;
		end
	end
clear f1 f2
end
y = sum(sum(error)) / (NT*NK);

