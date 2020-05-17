function WPrice = WeightedMethod(thet,S0,V0,K,S,V,T,A,invA,B)

% Heston Call price using the Weighted Method
% Requires a uniform grid for the stock price, volatility, and maturity
% INPUTS
%   thet = Weighing parameter
%   S0 = Spot price at which to price the call
%   V0 = variance at which to price the call
%   K  = Strike price
%   S = uniform grid for the stock price
%   V = uniform grid for the volatility
%   T = uniform grid for the maturity
%   A = "A" matrix for weighted method
%   invA = A inverse
%   B = "B" matrix for weighted method
% OUTPUT
%   y = 2-D interpolated European Call price

% Required vector lengths and time increment
NS = length(S);
NV = length(V);
NT = length(T);

% Initialize the u vector
u = zeros(NS*NV,1);

% U(0) vector - value of U(T) at maturity
Si = repmat(S',NV,1);
U = max(0, Si - K);

% Loop through the time increments, updating U(t) to U(t+1) at each step
for t=2:NT
	u = U;
	if thet==0
		U = B*u;         % Explicit Method
	elseif thet==1
		U = invA*u;      % Implicit Method
	else
		U = A\B*u;       % Weighted Method
	end
end

% Restack the U vector to output a matrix
UU = reshape(U,NS,NV);

% Interpolate to get the price at S0 and v0
WPrice = interp2(V,S,UU,V0,S0);
