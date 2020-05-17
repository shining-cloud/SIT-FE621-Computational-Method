function [x w] = GenerateGaussLobatto(N)
% Generate abscissas (x) and weights (w) for Gauss Lobatto integration

% The Legendre polynomial L and its derivative dL
n = N-1;
m = floor(n/2);
for k=0:m
 	L(k+1)  = (1/2^n)*(-1)^k*factorial(2*n-2*k)/factorial(k) ...
 		    / factorial(n-k)/factorial(n-2*k);
	dL(k+1) = (1/2^n)*(-1)^k*factorial(2*n-2*k)/factorial(k) ...
		    / factorial(n-k)/factorial(n-2*k)*(n-2*k);
end

% Fill in the blank powers of L and dL
% Legendre polynomial P and its derivative dP
for k=1:n+1
	if mod(k,2)==0
 		P(k)=0;
		dP(k) = dL(k/2);
	else
 		P(k) = L((k+1)/2);
		dP(k) = 0;
	end
end

% Abscissas are the roots of dP and augmented with the endpoints
x = sortrows(roots(dP));
x = [-1 x' 1]';

% Find the weights
w = zeros(n+1,1);
for j=2:n+1
	for k=1:n+1
		Poly(k) = P(k)*x(j)^(n+1-k);
	end
	w(j) = 2/N/(N-1)/sum(Poly)^2;
end
w(1) = 2/N/(N-1);
w(n+1) = w(1);

