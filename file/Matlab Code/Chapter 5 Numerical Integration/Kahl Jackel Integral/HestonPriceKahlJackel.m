function y = HestonPriceKahlJackel(PutCall,S,K,T,rf,q,param,trap,X,W)

% Kahl and Jackel (2005) Heston price, using an integration range of [0,1].
% From "Not-so-Complex Algorithms in the Heston model"

kappa = param(1);
theta = param(2);
sigma = param(3);
v0    = param(4);
rho   = param(5);
lambda = 0;


% Forward price
F = S*exp((rf-q)*T);

% Required quantities for integrand
Cinf = sqrt(1-rho^2)/sigma*(v0 + kappa*theta*T);

% Required quantities for f1
if (kappa - rho*sigma) ~= 0
    ImC1 = (exp((rho*sigma-kappa)*T)*theta*kappa + theta*kappa*((kappa-rho*sigma)*T-1))/2/(kappa-rho*sigma)^2;
    ImD1 = (1-exp((rho*sigma-kappa)*T))/2/(kappa-rho*sigma);
else
    ImC1 = kappa*theta*T^2/4;
    ImD1 = T/2;
end

% Required quantities for f2
ImC2 = -(exp(-kappa*T)*theta*kappa + theta*kappa*(kappa*T-1))/2/kappa^2;
ImD2 = - (1-exp(-kappa*T/2))/2/kappa;

for u = 1:length(X);
	% Transformation of the abscissa from [-1,1] to [0,1]
	x = 0.5*X(u) + 0.5;
	if x == 0;
		% Integrand at left abscissa 0
		y(u) = 0.5*(F-K);
	elseif x == 1 
		% Integrand at right abscissa 1
		f1 = log(F/K) + ImC1 + ImD1*v0;
		f2 = log(F/K) + ImC2 + ImD2*v0;
		y(u) = 0.5*(F-K) + (F*f1 - K*f1)/(pi*Cinf);
	else
		% Integrand at remaining abscissas
		f1 = HestonProb(-log(x)/Cinf,kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,1,trap);
		f2 = HestonProb(-log(x)/Cinf,kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,2,trap);
		y(u) = 0.5*(F-K) + (F*f1 - K*f2)/(x*pi*Cinf);
	end
	% Multiply by the weights
	z(u) = W(u)*y(u);
end

% The call price
KJCall = exp(-rf*T)*(1/2)*sum(z);

% The put price by put call parity
if strcmp(PutCall,'C')
	y = KJCall;
else
	y = KJCall - S*exp(-q*T) + K*exp(-rf*T);
end
