function [S V1 V2 SimPrice] = DHQuadExpSim(params,S0,Strike,Mat,r,q,T,N,PutCall)

% Euler, Milstein, Implicit Milstein, and Weighted Implicit Milstein
% discretization of the Double Heston model.
% INPUTS
%   scheme = Scheme for the variance and stock price
%            'E' Euler, 'M' Milstein, 'IM' Implicit Milstein, or
%            'WM' Weighted Explicit-Implicit Scheme
%   negvar = 'R' Reflection or 'T' Truncation of negative variances
%   params = Heston parameters
%   S0  = Spot Price
%   Mat = Maturity
%   r = riskf ree rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
%   alpha = Weight for explicit-implicit scheme
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances
%   F = Matrix of overriden negative variances

% Heston parameters
kappa1 = params(1);
theta1 = params(2);
sigma1 = params(3);
v01    = params(4);
rho1   = params(5);

kappa2 = params(6);
theta2 = params(7);
sigma2 = params(8);
v02    = params(9);
rho2   = params(10);

% Time increment
dt = Mat/T;

% Required quantities
K01 = -rho1*kappa1*theta1*dt/sigma1;
K11 = dt/2*(kappa1*rho1/sigma1 - 1/2) - rho1/sigma1;
K21 = dt/2*(kappa1*rho1/sigma1 - 1/2) + rho1/sigma1;
K31 = dt/2*(1-rho1^2);

K02 = -rho2*kappa2*theta2*dt/sigma2;
K12 = dt/2*(kappa2*rho2/sigma2 - 1/2) - rho2/sigma2;
K22 = dt/2*(kappa2*rho2/sigma2 - 1/2) + rho2/sigma2;
K32 = dt/2*(1-rho2^2);

% Initialize the volatilities and stock processes
V1 = zeros(T,N);
V2 = zeros(T,N);
S  = zeros(T,N);

% Starting values for the variance and stock processes
S(1,:)  = S0;          % Spot price 
V1(1,:) = v01;         % Heston initial variance
V2(1,:) = v02;         % Heston initial variance

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T;
		m1 = theta1 + (V1(t-1,i) - theta1)*exp(-kappa1*dt);
		s1 = V1(t-1,i)*sigma1^2*exp(-kappa1*dt)/kappa1*(1-exp(-kappa1*dt)) ...
            + theta1*sigma1^2/2/kappa1*(1-exp(-kappa1*dt))^2;
		phi1 = s1/(m1^2);
		p1 = (phi1-1)/(phi1+1);
		U1 = rand(1);
		if phi1 < 1/2
			b1 = sqrt(2/phi1 - 1 + sqrt(2/phi1*(2/phi1-1)));
			a1 = m1/(1+b1^2);
			Zv1 = normICDF(U1);
			V1(t,i) = a1*(b1+Zv1)^2;
		elseif phi1 >= 1/2
			if U1 <= p1
				V1(t,i) = 0;
			elseif U1 > p1
				beta1 = (1-p1)/m1;
				V1(t,i) = log((1-p1)/(1-U1))/beta1;
			end
		end

		m2 = theta2 + (V2(t-1,i) - theta2)*exp(-kappa2*dt);
		s2 = V2(t-1,i)*sigma2^2*exp(-kappa2*dt)/kappa2*(1-exp(-kappa2*dt)) ...
            + theta2*sigma2^2/2/kappa2*(1-exp(-kappa2*dt))^2;
		phi2 = s2/(m2^2);
		p2 = (phi2-1)/(phi2+1);
		U2 = rand(1);
		if phi2 < 1/2
			b2 = sqrt(2/phi2 - 1 + sqrt(2/phi2*(2/phi2-1)));
			a2 = m2/(1+b2^2);
			Zv2 = normICDF(U2);
			V2(t,i) = a2*(b2+Zv2)^2;
		elseif phi2 >= 1/2
			if U2 <= p2
				V2(t,i) = 0;
			elseif U2 > p2
				beta2 = (1-p2)/m2;
				V2(t,i) = log((1-p2)/(1-U2))/beta2;
			end
		end
		% Predictor-Corrector for the stock price
		B1 = randn(1);
		B2 = randn(1);
		logS = log(exp(-r*t*dt)*S(t-1,i)) ...
			 + K01 + K11*V1(t-1,i) + K21*V1(t,i) + sqrt(K31*(V1(t,i)+V1(t-1,i)))*B1 ...
			 + K02 + K12*V2(t-1,i) + K22*V2(t,i) + sqrt(K32*(V2(t,i)+V2(t-1,i)))*B2;
		S(t,i) = exp(logS)*exp(r*(t+1)*dt);
	end
end

% Terminal stock prices
ST = S(end,:);

% Payoff vectors
if strcmp(PutCall,'C')
	Payoff = max(ST - Strike,0);
elseif strcmp(PutCall,'P')
	Payoff = max(Strike - ST,0);
end
	
% Simulated price
SimPrice = exp(-r*Mat)*mean(Payoff);
