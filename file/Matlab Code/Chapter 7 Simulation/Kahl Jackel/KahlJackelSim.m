function [S V F] = KahlJackelSim(scheme,negvar,params,S0,Mat,r,q,T,N)

% Kahl-Jackel discretization of the Heston model.
% Runs the IJK scheme, Pathwise scheme, Balanced Implicit Scheme
% INPUTS
%   scheme = 'IJK' for IJK, 'PW' for pathwise, 'B' Balanced 
%   negvar = 'R' Reflection or 'T' Truncation of negative variances
%   params = Heston parameters
%   S0  = Spot Price
%   Mat = Maturity
%   r = riskf ree rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances
%   F = Matrix of overriden negative variances

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

% Flags for negative variances
F = 0;

% Starting values for the variance and stock processes
S(1,:) = S0;       % Spot price 
V(1,:) = v0;       % Heston v0 initial variance 

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T;
		% Generate two dependent N(0,1) variables with correlation rho
		Zv = randn(1);
		Zs = rho*Zv + sqrt(1-rho^2)*randn(1);
		if strcmp(scheme,'IJK')
			% Implicit Milstein for the variance.
			V(t,i) = (V(t-1,i) + kappa*theta*dt ...
				+ sigma*sqrt(V(t-1,i)*dt)*Zv...
				+ sigma^2*dt*(Zv^2-1)/4) / (1+kappa*dt);

			% Apply the full truncation or reflection scheme to the variance
			if V(t,i) <= 0
				F = F+1;
				if strcmp(negvar,'R')          % Reflection: take -V
					V(t,i) = abs(V(t,i));
				elseif strcmp(negvar,'T')
					V(t,i) = max(0, V(t,i));   % Truncation: take max(0,V)
				end
			end

			% IJK discretization scheme for the log stock prices
			S(t,i) = S(t-1,i)*exp( ...
				+ (r-q-(V(t,i)+V(t-1,i))/4)*dt...
				+ rho*sqrt(V(t-1,i)*dt)*Zv ...
				+ (1/2)*(sqrt(V(t-1,i))+sqrt(V(t,i)))*(Zs-rho*Zv)*sqrt(dt)...
				+ rho*sigma*dt*(Zv^2-1)/2);

		elseif strcmp(scheme,'PW')
			% Pathwise Adapated Linearization Quadratic for the variance.
			% Note: Bn = dZ/dt = Z*sqrt(dt)/dt = Z/sqrt(dt);
			theta2 = theta - sigma^2/4/kappa;
			Bn = Zv/sqrt(dt);
			V(t,i) = V(t-1,i)...
				+ (kappa*(theta2-V(t-1,i)) + sigma*Bn*sqrt(V(t-1,i)))*dt...
				*(1 + (sigma*Bn-2*kappa*sqrt(V(t-1,i)))*dt/4/sqrt(V(t-1,i)));

			% Apply the full truncation or reflection scheme to the variance
			if V(t,i) <= 0
				F = F+1;
				if strcmp(negvar,'R')          % Reflection: take -V
					V(t,i) = abs(V(t,i));
				elseif strcmp(negvar,'T')
					V(t,i) = max(0, V(t,i));   % Truncation: take max(0,V)
				end
			end

			% Euler/Milstein discretization scheme for the log stock prices
			S(t,i) = S(t-1,i)*exp((r-q-V(t-1,i)/2)*dt + sqrt(V(t-1,i)*dt)*Zs);
            
		elseif strcmp(scheme,'B')
			% Balanced Implicit scheme for the variance
%			absdW = sqrt(dt)*abs(Zv);
			C = kappa*dt + sigma/sqrt(V(t-1,i))*sqrt(dt)*abs(Zv);
			V(t,i) = (V(t-1,i)*(1+C) + kappa*(theta-V(t-1,i))*dt ...
				   + sigma*sqrt(V(t-1,i)*dt)*Zv) / (1+C);
               
			% Euler/Milstein discretization scheme for the log stock prices
			S(t,i) = S(t-1,i)*exp((r-q-V(t-1,i)/2)*dt + sqrt(V(t-1,i)*dt)*Zs);
		end
	end
end


