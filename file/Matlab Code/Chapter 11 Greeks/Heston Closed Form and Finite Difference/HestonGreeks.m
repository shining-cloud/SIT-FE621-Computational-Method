function y = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,Greek)

% Heston (1993) Greeks by Gauss-Laguerre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T = Time to maturity.
%   r = Risk free rate.
%   kappa  = Heston parameter: mean reversion speed.
%   theta  = Heston parameter: mean reversion level.
%   sigma  = Heston parameter: volatility of vol
%   lambda = Heston parameter: risk.
%   v0     = Heston parameter: initial variance.
%   rho    = Heston parameter: correlation
%   trap:  1 = "Little Trap" formulation
%          0 = Original Heston formulation
%   x = Gauss Laguerre abscissas
%   w = Gauss Laguerre weights
% Greek = 'Delta', 'Rho', 'Gamma', 'Theta', 'Vega1'
%         'Vega2', 'Vanna', 'Volga', 'Kappa','Sigma','Corr'
% OUTPUT -------------------------------------------------------
%   The Heston Greek

if strcmp(Greek,'Delta')
	for k=1:length(x);
		% P1
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Delta');
	end
	P1 = 1/2 + 1/pi*sum(int1);
	if strcmp(PutCall,'C')
		y = exp(-q*T)*P1;
	else
		y = exp(-q*T)*(P1 - 1);
	end
elseif strcmp(Greek,'Gamma')
	for k=1:length(x)
		% der(P1^2)/der(S^2)
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Gamma');
	end
	y = exp(-q*T)/pi/S * sum(int1);
elseif strcmp(Greek,'Rho')
	for k=1:length(x)
		int2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,'Rho');
	end
	P2 = 1/2 + 1/pi*sum(int2);
	if strcmp(PutCall,'C')
		y = K*exp(-r*T)*T*P2;
	else
		y = K*exp(-r*T)*T*(P2-1);
	end
elseif strcmp(Greek,'Theta')
	for k=1:length(x)
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Delta');
		int2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,'Delta');
		der1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Theta');
		der2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,'Theta');
	end
	P1 = 1/2 + 1/pi*sum(int1);
	P2 = 1/2 + 1/pi*sum(int2);
	dP1 = 1/pi*sum(der1);
	dP2 = 1/pi*sum(der2);
	theta = -S*exp(-q*T)*(-q*P1 + dP1) + K*exp(-r*T)*(dP2-r*P2);
	if strcmp(PutCall,'C')
		y = theta;
	elseif strcmp(PutCall,'P')
		y = theta + K*r*exp(-r*T) - q*S*exp(-q*T);
	end
elseif strcmp(Greek,'Vega1')
	for k=1:length(x)
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Vega1');
		int2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,'Vega1');
	end
	dP1 = 1/pi*sum(int1)*2*sqrt(v0);
	dP2 = 1/pi*sum(int2)*2*sqrt(v0);
	y = S*exp(-q*T)*dP1 - K*exp(-r*T)*dP2;
elseif strcmp(Greek,'Vega2')
	for k=1:length(x)
 		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Vega2');
 		int2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,'Vega2');
	end
	dP1 = 1/pi*sum(int1)*2*sqrt(theta);
	dP2 = 1/pi*sum(int2)*2*sqrt(theta);
	y = S*exp(-q*T)*dP1 - K*exp(-r*T)*dP2;
elseif strcmp(Greek,'Vanna');
	for k=1:length(x);
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Vega1');
	end
	dP1 = 1/pi*sum(int1)*2*sqrt(v0);
	y = exp(-q*T)*dP1;
elseif strcmp(Greek,'Volga')
	% Calculate Vega1
	for k=1:length(x)
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Vega1');
		int2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,'Vega1');
	end
	dP1 = 1/pi*sum(int1)*2*sqrt(v0);
	dP2 = 1/pi*sum(int2)*2*sqrt(v0);
	Vega1 = S*exp(-q*T)*dP1 - K*exp(-r*T)*dP2;
	% Calculate second-order derivative of the call price w.r.t. v0
	for k=1:length(x);
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,'Volga');
		int2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,'Volga');
	end
	dP1 = 1/pi*sum(int1);
	dP2 = 1/pi*sum(int2);
	dC2 = S*exp(-q*T)*dP1 - K*exp(-r*T)*dP2;
	% Calculate Volga
	y = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1/4/v0);
elseif strcmp(Greek,'Corr') | strcmp(Greek,'Sigma') | strcmp(Greek,'Kappa');
	for k=1:length(x)
		int1(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap,Greek);
		int2(k) = w(k)*HestonGreeksProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap,Greek);
	end
	dP1 = 1/pi*sum(int1);
	dP2 = 1/pi*sum(int2);
	y = S*exp(-q*T)*dP1 - K*exp(-r*T)*dP2;
end


	