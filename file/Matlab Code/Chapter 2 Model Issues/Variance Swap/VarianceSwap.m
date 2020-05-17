function y = VarianceSwap(KCI,CallVI,KPI,PutVI,S,T,rf,q)

% Demeterfi et al (1999) fair strike of a variance swap
% INPUTS
%   KCI    = Grid of call strikes
%   CallVI = Interpolated call implied vol along the grid KCI
%   KPI    = Grid of put strikes
%   PutVI  = Interpolated put implied vol along the grid CallKI
%   S = Spot price
%   T = Maturity
%   rf = interest rate
%   q = dividend yield

% Do the required calculations on calls.  
n = length(CallVI);

% Take ATM as the boundary point
Sb = S;

% Auxiliary Functions
f = @(S, Sb, T) 2/T*((S - Sb) / Sb - log(S/Sb));
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));

% Start the replication algorithm
for i=1:n-1
	Temp(i) = (f(KCI(i+1),Sb,T) - f(KCI(i),Sb,T)) / (KCI(i+1) - KCI(i));
	if i==1 
		CallWeight(1) = Temp(1);
	end
	CallValue(i) = BSC(S,KCI(i),rf,0,CallVI(i),T);
	if i>1
		CallWeight(i) = Temp(i) - Temp(i-1);
	end;
    CallContrib(i) = CallValue(i)*CallWeight(i);
end
Pi1 = sum(CallContrib);

% Do the calculations on puts. Flip the Vectors for Convenience
n = length(PutVI);
KPI   = fliplr(KPI);
PutVI = fliplr(PutVI);

for i=1:n-1
    Temp2(i) = (f(KPI(i+1),Sb,T) - f(KPI(i),Sb,T)) / (KPI(i) - KPI(i+1));
	if i==1
		PutWeight(1) = Temp2(1);
	end;
	PutValue(i) = BSP(S,KPI(i),rf,0,PutVI(i),T);
	if i>1
		PutWeight(i) = Temp2(i) - Temp2(i-1);
	end
    PutContrib(i) = PutValue(i) * PutWeight(i);
end
Pi2 = sum(PutContrib);

% Total cost of the portfolio
Pi_CP = Pi1 + Pi2;


% Results of the replication
% Estimate of fair variance
Kvar = 2/T*(rf*T - (S/Sb*exp(rf*T) - 1) - log(Sb/S)) + exp(rf*T)*Pi_CP;

y = Kvar;
