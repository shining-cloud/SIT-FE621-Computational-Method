function y = TrinomialTree(Spot,K,r,q,v,T,N,PutCall,EuroAmer)

% Trinomial tree parameters and probabilities.
dt = T/N;
u = exp(v*sqrt(2*dt));
d = 1/u;
pu = (exp((r-q)*dt/2) - exp(-v*sqrt(dt/2)))^2/(exp(v*sqrt(dt/2)) - exp(-v*sqrt(dt/2)))^2;
pd = (exp(v*sqrt(dt/2)) - exp((r-q)*dt/2))^2/(exp(v*sqrt(dt/2)) - exp(-v*sqrt(dt/2)))^2;
pm = 1 - pu - pd;

% Initialize the stock prices
S = zeros(2*N+1,N+1);
S(N+1,1) = Spot;

% Calculate all the stock prices.
for j=2:N+1
    for i=N-j+2:N+j
        S(i,j) = Spot*u^(N+1-i);
	end
end

% Initialize the option prices.
V = zeros(2*N+1,N+1);

% Calculate terminal option prices.
switch PutCall
	case 'C'
		V(:,N+1) = max(S(:,N+1) - K, 0);
	case 'P'
		V(:,N+1) = max(K - S(:,N+1), 0);
end

% Calculate remaining entries for Calls and Puts
for j=N:-1:1
	for i=N-j+2:N+j
		switch EuroAmer
			case 'A'
				if strcmp(PutCall, 'C')
					V(i,j) = max(S(i,j) - K, exp(-r*dt)*(pu*V(i-1,j+1) + pm*V(i,j+1) + pd*V(i+1,j+1)));
				else
					V(i,j) = max(K - S(i,j), exp(-r*dt)*(pu*V(i-1,j+1) + pm*V(i,j+1) + pd*V(i+1,j+1)));
				end
			case 'E'
				V(i,j) = exp(-r*dt)*(pu*V(i-1,j+1) + pm*V(i,j+1) + pd*V(i+1,j+1));
		end
	end
end

% Option price is at the first node.
y = V(N+1,1);

	