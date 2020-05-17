function [x w] = GenerateGaussLegendre(n)
% Generate abscissas (x) and weights (w) for Gauss Legendre integration

% The Legendre polynomial
m = floor(n/2);
for k=0:m
	L(k+1) = (1/2^n)*(-1)^k*factorial(2*n-2*k)/factorial(k)/factorial(n-k)/factorial(n-2*k);
end

% Fill in the blank powers of L
for k=1:n+1
	if mod(k,2)==0
		P(k)=0;
	else
		P(k) = L((k+1)/2);
	end
end

%Find the roots
x = sortrows(roots(P));

% Find the weights
w = zeros(n,1);
for j=1:n;
	% The coefficients of the derivative of the Legrenge polynomial
	for k=0:m
		dC(k+1,j) = (1/2^n)*(-1)^k*factorial(2*n-2*k)/factorial(k)/factorial(n-k)/factorial(n-2*k)*(n-2*k)*x(j)^(n-2*k-1);
	end
	% The weight w(j)
	w(j) = 2/(1-x(j)^2)/sum(dC(:,j))^2;
end
