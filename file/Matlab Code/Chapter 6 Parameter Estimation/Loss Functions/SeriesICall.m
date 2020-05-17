function y = SeriesICall(S,K,rf,q,T,v0,rho,theta,kappa,sigma)

% Black Scholes call price and Black Scholes vega
% Note that v = variance, not volatility
BSC = @(S,K,rf,q,v,T) (S*exp(-q*T)*normcdf((log(S/K) + (rf-q+v/2)*T)/sqrt(v*T)) - K*exp(-rf*T)*normcdf((log(S/K) - (rf-q+v/2)*T)/sqrt(v*T)));
BSV = @(S,K,rf,q,v,T) (sqrt(T/8/pi/v)*S*exp(-q*T)*exp(-0.5*((log(S/K) + (rf-q+v/2)*T)/sqrt(v*T))^2));

% The "J" integrals
J1 = J(rho,theta,kappa,T,v0,1);
J3 = J(rho,theta,kappa,T,v0,3);
J4 = J(rho,theta,kappa,T,v0,4);

% Time average of the deterministic volatility
v = theta + (v0-theta)*(1-exp(-kappa*T))/(kappa*T);

% X and Z required for the ratios of Black Scholes derivatives
X = log(S*exp((rf-q)*T)/K);
Z = v*T;

% The ratios of Black Scholes derivatives
R20 = R(2,0,X,Z,T);
R11 = R(1,1,X,Z,T);
R12 = R(1,2,X,Z,T);
R22 = R(2,2,X,Z,T);

% Black Scholes call price evaluated at v
c = BSC(S,K,rf,q,v,T);

% Black Scholes vega evaluated at v
cv = BSV(S,K,rf,q,v,T);

% Series I volatility of volatility expansion for the implied volatility
iv = c ...
   + sigma/T*J1*R11*cv ...
   + sigma^2*(J3*R20/T^2 + J4*R12/T + J1^2*R22/T^2/2)*cv;

% Error control on the volatility
if iv <= 0
	iv = 0.01;
end

% Heston Price as Black-Scholes, with the Lewis IV as the volatility parameter
y = BSC(S,K,rf,q,iv,T);

		 