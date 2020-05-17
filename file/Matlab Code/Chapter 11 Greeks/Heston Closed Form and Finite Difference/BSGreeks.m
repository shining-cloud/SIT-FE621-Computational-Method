function y = BSGreeks(PutCall,S,K,r,q,T,sigma,Greek)

% Calculates the Black-Scholes Greeks

d1 = (log(S/K) + (r-q+sigma^2/2)*T)/sigma/sqrt(T);
d2 = d1 - sigma*sqrt(T);
Nd1 = normcdf(d1);
Nd2 = normcdf(d2);
nd1 = normpdf(d1);

if strcmp(Greek,'Delta')
	if strcmp(PutCall,'C')
		y = exp(-q*T)*Nd1;
	else
		y = exp(-q*T)*(Nd1 - 1);
	end
elseif strcmp(Greek,'Gamma')
	y = nd1*exp(-q*T)/S/sigma/sqrt(T);
elseif strcmp(Greek,'Vega')
	y = S*exp(-q*T)*sqrt(T)*nd1;
elseif strcmp(Greek,'Rho')
	if strcmp(PutCall,'C')
		y = K*T*exp(-r*T)*Nd2;
	else
		y = K*T*exp(-r*T)*(Nd2-1);
	end
elseif strcmp(Greek,'Theta')
	if strcmp(PutCall,'C')
		y = -exp(-q*T)*S*nd1*sigma/2/sqrt(T) - r*K*exp(-r*T)*Nd2 + q*S*exp(-q*T)*Nd1;
	else
		y = -exp(-q*T)*S*nd1*sigma/2/sqrt(T) + r*K*exp(-r*T)*(1-Nd2) - q*S*exp(-q*T)*(1-Nd1);
	end
elseif strcmp(Greek,'Vanna')
	y = -exp(-q*T)*nd1*d2/sigma;
elseif strcmp(Greek,'Volga')
	y = S*exp(-q*T)*nd1*sqrt(T)*d1*d2/sigma;
end

	