function [U u] = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer)

% Finite differences for the Heston PDE for a European Call
% Uses uneven grid sizes
% In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
% in the Heston Modelo with Correlation" Int J of Num Analysis and Modeling, 2010.
% INPUTS
%    params = 6x1 vector of Heston parameters
%    K = Strike price
%    r = risk free rate
%    q = dividend yield
%    S = vector for stock price grid
%    V = vector for volatility grid
%    T = vector for maturity grid
% OUTPUTS
%    U = U(S,v) 2-D array of size (nS+1)x(nV+1)x(nT+1) for the call price

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);

% Grid measurements
NS = length(S);
NV = length(V);
NT = length(T);
Smin = S(1); Smax = S(NS);
Vmin = V(1); Vmax = V(NV);
Tmin = T(1); Tmax = T(NT);
dt = (Tmax-Tmin)/(NT-1);

% Initialize the 2-D grid with zeros
U = zeros(NS,NV);

% Temporary grid for previous time steps
u = zeros(NS,NV);

%% Solve the PDE
% Round each value of U(S,v) at each step
% Boundary condition for t=maturity
for s=1:NS
    if strcmp(PutCall,'C')
        U(s,:) = max(S(s) - K, 0);
    elseif strcmp(PutCall,'P')
        U(s,:) = max(K - S(s), 0);
    end
end

% Go through the times
for t=1:NT-1
    % Boundary condition for Smin and Smax
    U(1,:) = 0;
    if strcmp(PutCall,'C')
        U(NS,:) = max(0, Smax - K); % Galiotos uses U(NS-1,v) + ds;
        U(:,NV) = max(0, S - K); % Galiotos uses U(s,NV-1);
    elseif strcmp(PutCall,'P')
        U(NS,:) = max(0, K - Smax); % Galiotos uses U(NS-1,v) + ds;
        U(:,NV) = max(0, K - S); % Galiotos uses U(s,NV-1);
    end
    % Update the temporary grid u(s,t) with the boundary conditions
    u = U;
    % Boundary condition for Vmin.
    % Previous time step values are in the temporary grid u(s,t)
    for s=2:NS-1
        derV = (u(s,2)   - u(s,1))   / (V(2)-V(1));          % Forward difference
        derS = (u(s+1,1) - u(s-1,1)) / (S(s+1)-S(s-1));      % Central difference
        LHS = - r*u(s,1) + (r-q)*S(s)*derS + kappa*theta*derV;
        U(s,1) = LHS*dt + u(s,1);
    end
    % Update the temporary grid u(s,t) with the boundary conditions
    u = U;
    % Interior points of the grid (non boundary).
    % Previous time step values are in the temporary grid u(s,t)
    for s=2:NS-1
        for v=2:NV-1
            derS  =  (u(s+1,v) - u(s-1,v)) / (S(s+1)-S(s-1));  % Central difference for dU/dS
            derV  =  (u(s,v+1) - u(s,v-1)) / (V(v+1)-V(v-1));  % Central difference for dU/dV
            derSS = ((u(s+1,v) - u(s,v))   / (S(s+1)-S(s)) - (u(s,v) - u(s-1,v))/(S(s)-S(s-1))) / (S(s+1)-S(s));  % d2U/dS2
            derVV = ((u(s,v+1) - u(s,v))   / (V(v+1)-V(v)) - (u(s,v) - u(s,v-1))/(V(v)-V(v-1))) / (V(v+1)-V(v));  % d2U/dV2
            derSV =  (u(s+1,v+1) - u(s-1,v+1) - U(s+1,v-1) + U(s-1,v-1)) / (S(s+1)-S(s-1)) / (V(v+1)-V(v-1));     % d2U/dSdV
            L = 0.5*V(v)*S(s)^2*derSS + rho*sigma*V(v)*S(s)*derSV...
                + 0.5*sigma^2*V(v)*derVV - r*u(s,v)...
                + (r-q)*S(s)*derS + kappa*(theta-V(v))*derV;
            U(s,v) = L*dt + u(s,v);
        end
    end
    if strcmp(EuroAmer,'A')
        for s=1:NS
            if strcmp(PutCall,'C')
                U(s,:) = max(U(s,:), S(s) - K);
            elseif strcmp(PutCall,'P')
                U(s,:) = max(U(s,:), K - S(s));
            end
        end
    end
end


