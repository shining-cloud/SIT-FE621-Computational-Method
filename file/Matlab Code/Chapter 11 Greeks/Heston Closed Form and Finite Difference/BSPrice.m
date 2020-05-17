function y = BSPrice(PutCall,S,K,r,q,T,sigma);

d1 = (log(S/K) + (r-q+sigma^2/2)*T)/sigma/sqrt(T);
d2 = d1 - sigma*sqrt(T);
Nd1 = normcdf(d1);
Nd2 = normcdf(d2);

BSCall = exp(-q*T)*S*Nd1 - exp(-r*T)*K*Nd2;

if strcmp(PutCall,'C')
	y = BSCall;
else
	y = BSCall + K*exp(-r*T) - S*exp(-q*T);
end

	
	
	
	