function [S V] = QESim(params,gamma1,gamma2,S0,Mat,r,q,T,N,MC,phic,icdf)

% Quadratic Exponential (QE) Algorithm
% With martingale correction on the stock price simulation.
% INPUTS
%   params = Heston parameters
%   gamma1 = first gamma parameter
%   gamma2 = second gamma parameter
%   S0  = Spot Price
%   Mat = Maturity
%   r = riskf ree rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
%   MC = 1:Martingale correction, 0:No correction
%   phic = cut-off for low and high regimes for v(t)
%   icdf = 1=Matlab algorithm
%          2=Wichura algorithm
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);

% Time increment
dt = Mat/T;

% Initialize the variance and stock processes
V = zeros(T,N);
S = zeros(T,N);

% Starting values for the variance and stock processes
S(1,:) = S0;       % Spot price 
V(1,:) = v0;       % Heston v0 initial variance 

% Parameters for the QE scheme
E = exp(-kappa*dt);  % Simplify the coding
K0 = -kappa*rho*theta*dt/sigma;
K1 = (kappa*rho/sigma - 1/2)*gamma1*dt - rho/sigma;
K2 = (kappa*rho/sigma - 1/2)*gamma2*dt + rho/sigma;
K3 = gamma1*dt*(1-rho^2);
K4 = gamma2*dt*(1-rho^2);
A = K2 + K4/2;

% Loop through the simulation runs
for i=1:N;
	% Loop through the time increments
	for t=2:T;
		% Generate two dependent N(0,1) random variables
		Zv = randn(1);
		Zs = rho*Zv + sqrt(1-rho^2)*randn(1); 

		% QE Agorithm
		m = theta + (V(t-1,i) - theta)*E;
		s2 = V(t-1,i)*sigma^2*E/kappa*(1-E) + theta*sigma^2/2/kappa*(1-E)^2;
		phi = s2/m^2;
		Uv = rand;
		if phi <= phic
			b = sqrt(2/phi - 1 + sqrt(2/phi*(2/phi-1)));
			a = m/(1+b^2);
			if icdf==1
				Zv = norminv(Uv);
			elseif icdf==2
				Zv = normICDF(Uv);
			end
			V(t,i) = a*(b + Zv)^2;
			% Martingale correction: Define new K0 if possible
			if (MC==1) & A<(1/(2*a));
                M = exp(A*b^2*a/(1-2*A*a))/sqrt(1-2*A*a);
                K0 = -log(M) - (K1+0.5*K3)*V(t,i);
			end
			S(t,i) = S(t-1,i)*exp((r-q)*dt + K0 + K1*V(t-1,i) + K2*V(t,i) + ...
				                  sqrt(K3*V(t-1,i) + K4*V(t,i))*Zs);
		else
			p = (phi-1)/(phi+1);
			beta = (1-p)/m;
			if (0<=Uv) & (Uv<=p);
				phiinv = 0;
			elseif (p<Uv) & (Uv<=1);
				phiinv = 1/beta*log((1-p)/(1-Uv));
			end
			V(t,i) = phiinv;
			% Martingale correction: Define new K0 if possible
			if MC==1 & A<beta;
                M = p + beta*(1-p)/(beta-A);
                K0 = -log(M) - (K1+0.5*K3)*V(t,i);
			end			
			S(t,i) = S(t-1,i)*exp((r-q)*dt + K0 + K1*V(t-1,i) + K2*V(t,i) + ...
				                  sqrt(K3*V(t-1,i) + K4*V(t,i))*Zs);			
		end;
	end
end

