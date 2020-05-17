function U = HestonExplicitPDE(params,K,r,q,S,V,T)

% Finite differences for the Heston PDE for a European Call
% Uses even grid sizes
% In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
% in the Heston Modelo with Correlation" Int J of Num Analysis and Modeling, 2010.
% Thesis by Sensi Li and paper by Vassilis Galiotos
% INPUTS
%    params = 6x1 vector of Heston parameters
%    K = Strike price
%    r = risk free rate
%    q = Dividend yield
%    S = vector for stock price grid
%    V = vector for volatility grid
%    T = vector for maturity grid
% OUTPUT
%    U = U(S,v) 2-D array of size (nS+1)x(nV+1) for the call price

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);

% Length of stock price, volatility, and maturity
NS = length(S);
NV = length(V);
NT = length(T);
Smin = S(1);  Smax = S(NS);
Vmin = V(1);  Vmax = V(NV);
Tmin = T(1);  Tmax = T(NT);

% Increment for Stock Price, Volatility, and Maturity
ds = (Smax-Smin)/(NS-1);
dv = (Vmax-Vmin)/(NV-1);
dt = (Tmax-Tmin)/(NT-1);

% Initialize the 2-D grid with zeros
U = zeros(NS,NV);

% Temporary grid for previous time steps
u = zeros(NS,NV);

% Solve the PDE
% Round each value of U(S,v,t) at each step
% Boundary condition for t = Maturity
for s=1:NS
	for v=1:NV
		U(s,v) = max(S(s) - K, 0);
	end
end

for t=1:NT-1
	% Boundary condition for Smin and Smax
	for v=1:NV-1
		U(1,v) = 0;
		U(NS,v) = max(0, Smax - K); % Galiotos uses U(NS-1,v) + ds;
	end
	% Boundary condition for Vmax
	for s=1:NS
		U(s,NV) = max(0, S(s) - K); % Galiotos uses U(s,NV-1);
	end
	% Update the temporary grid u(s,t) with the boundary conditions
	u = U;
	% Boundary condition for Vmin.
	% Previous time step values are in the temporary grid u(s,t)
	for s=2:NS-1
		DerV = (u(s,2) - u(s,1)) / dv;                	% PDE Points on the middle of the grid (non boundary)
		DerS = (u(s+1,1) - u(s-1,1))/2/ds;              % Central difference for dU/dS
		U(s,1) = u(s,1)*(1 - r*dt - kappa*theta*dt/dv)...
			   + dt*(0.5*(r-q)*s*u(s+1,1) - u(s-1,1)) ...
			   + kappa*theta*dt/dv*u(s,2);
	end
	% Update the temporary grid u(s,t) with the boundary conditions
	u = U;
	% Interior points of the grid (non boundary).
	% Previous time step values are in the temporary grid u(s,t)
	for s=2:NS-1
		for v=2:NV-1
			A = (1 - dt*(s-1)^2*(v-1)*dv - sigma^2*(v-1)*dt/dv - r*dt);
			B = (1/2*dt*(s-1)^2*(v-1)*dv - 1/2*dt*(r-q)*(s-1));
			C = (1/2*dt*(s-1)^2*(v-1)*dv + 1/2*dt*(r-q)*(s-1));
			D = (1/2*dt*sigma^2*(v-1)/dv - 1/2*dt*kappa*(theta-(v-1)*dv)/dv);
			E = (1/2*dt*sigma^2*(v-1)/dv + 1/2*dt*kappa*(theta-(v-1)*dv)/dv);
			F = 1/4*dt*sigma*(s-1)*(v-1);
			U(s,v) = A*u(s,v) + B*u(s-1,v) + C*u(s+1,v)...   % The PDE
				   + D*u(s,v-1) + E*u(s,v+1)...
				   + F*(u(s+1,v+1)+u(s-1,v-1)-u(s-1,v+1)-u(s+1,v-1));
		end
	end
end

