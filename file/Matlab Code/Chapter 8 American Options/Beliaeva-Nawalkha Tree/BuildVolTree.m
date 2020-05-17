function [X V RBound M] = BuildVolTree(kappa,theta,sigma,V0,dt,NT,threshold)

% Creates the Beliaeva-Nawalkha tree for the variance process
% INPUTS
%   kappa = Heston parameter for mean reversion speed
%   theta = Heston parameter for mean reversion level
%   sigma = Heston parameter for vol of vol
%   V0     = Heston parameter for initial variance
%   rf = Risk free rate
%   dt = time increment
%   NT = Number of time steps
%   threshold = threshold for which V can be zero
% OUTPUTS
%   X  = the tree for transformed variance
%   V  = the tree for the original variance
%   RBound = the row where the V(n,t) variances are zero
%   M = row where volatility starts at time 1

%% Initial quantities
X0 = 2*sqrt(V0)/sigma;

% Equations (24) through (26)
be = X0/sqrt(dt)/floor(X0/sqrt(1.5*dt));
bc = X0/sqrt(dt)/floor(X0/sqrt(1.5*dt)+1);
if abs(bc-sqrt(1.5)) < abs(be-sqrt(1.5))
	b = bc;
else
	b = be;
end

if (b<1) || (b>sqrt(2))
    fprintf('Warning b = %5.3f but should be 1 < b < 1.4142\n',b);
end

% Initialize the X-matrix and V-matrix
NR = 2*NT-1;
X = zeros(NR,NT);
V = zeros(NR,NT);
M = (NR+1)/2;
X(M,1) = X0;


%% Time 1 node for X
% Equations (22) and (30) 
muX = 1/X(M,1)*(0.5*kappa*(4*theta/sigma^2-X(M,1)^2)-0.5);
J = floor(muX*sqrt(dt)/b + 1/b^2);

% Initialize the first nodes
X(M-1,2) = X(M,1) + b*(J+1)*sqrt(dt);
X(M+0,2) = X(M,1) + b*(J+0)*sqrt(dt);
X(M+1,2) = X(M,1) + b*(J-1)*sqrt(dt);


%% Remaining nodes for X
for t=2:M-1
    for n=M-t+1:M+t-1
        muX = 1/X(n,t)*(0.5*kappa*(4*theta/sigma^2-X(n,t)^2)-0.5);
        J = floor(muX*sqrt(dt)/b + 1/b^2);
        if X(n,t)>threshold && X(n,t)^2*sigma^2/4>threshold
            % Case 1: Nodes where X > 0 -- Equation (27)
            X(n-1,t+1) = X(n,t) + b*(J+1)*sqrt(dt);       % up
            X(n+0,t+1) = X(n,t) + b*(J+0)*sqrt(dt);       % middle
            X(n+1,t+1) = X(n,t) + b*(J-1)*sqrt(dt);       % down
            % Case 2: Nodes where X = 0
        else
              X(n-1,t+1) = X(n-1,t);
              X(n-0,t+1) = X(n-0,t);
              X(n+1,t+1) = X(n+1,t);
%               X(n-1,t+1) = X(n,t) + b*(J+1)*sqrt(dt);
%               X(n+0,t+1) = X(n+0,t);
%               X(n+1,t+1) = X(n+1,t);
        end
    end
end


%% Identify row of the tree where X = 0
RBound = 1;
while X(RBound,M)>=threshold && RBound < NR;
    RBound = RBound+1;
end


%% Build the volatility tree V(n,t)
for t=1:NT
    for n=1:NR
        V(n,t) = X(n,t)^2*sigma^2/4;
        if V(n,t) < threshold
            V(n,t) = 0;
        end
        if X(n,t) < threshold
            X(n,t) = 0;
        end
    end
end

