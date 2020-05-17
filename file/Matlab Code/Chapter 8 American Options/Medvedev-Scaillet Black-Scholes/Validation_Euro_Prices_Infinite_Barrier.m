% Verifies that choosing y = inf yields European prices

clc; clear;

% Black Scholes European price
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));

% Option settings
S = 40;
r = 0.0488;
q = 0.00;

% Trinomial tree settings
N = 500;
PutCall = 'P';
EuroAmer = 'A';

% Table 1 settings
sigma = [0.2 0.2 0.2 0.3 0.3 0.3 0.4 0.4 0.4];
T = repmat([1/12 1/3 7/12],1,3);

EuroPut35 = zeros(3,9);     % K = 35
EuroPut40 = zeros(3,9);     % K = 40
EuroPut45 = zeros(3,9);     % K = 45

% Set y = inf to obtain European prices
y = 1e50;

%% Find the tree and M-S prices for K = 35
K = 35;
for i=1:9
    EuroPut35(1,i) = BSP(S,K,r,q,sigma(i),T(i));
    theta = log(K/S)/sigma(i)/sqrt(T(i));
    EuroPut35(2,i) = MSPutBS(y,theta,K,sigma(i),r,q,T(i));
end

%% Find the tree and M-S prices for K = 40
K = 40;
for i=1:9
    EuroPut40(1,i) = BSP(S,K,r,q,sigma(i),T(i));
    theta = log(K/S)/sigma(i)/sqrt(T(i));
    EuroPut40(2,i) = MSPutBS(y,theta,K,sigma(i),r,q,T(i));
end

%% Find the tree and M-S prices for K = 45
K = 45;
for i=1:9
    EuroPut45(1,i) = BSP(S,K,r,q,sigma(i),T(i));
    theta = log(K/S)/sigma(i)/sqrt(T(i));
    EuroPut45(2,i) = MSPutBS(y,theta,K,sigma(i),r,q,T(i));
end


%% Output the results
fmtstring = repmat('%8.5f',1,9);
fprintf('                sigma = 0.2             sigma = 0.3             sigma = 0.4\n')
fprintf('           ----------------------  ----------------------  ------------------------\n');
fprintf('K = 35     T=1/12  T=1/3   T=7/12  T=1/12  T=1/3   T=7/12  T=1/12  T=1/3   T=7/12  \n');
fprintf('-----------------------------------------------------------------------------------\n')
fprintf(['ClosedPut' fmtstring '\n'], EuroPut35(1,:))
fprintf(['MSPut    ' fmtstring '\n'], EuroPut35(2,:))
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('K = 40\n')
fprintf('-----------------------------------------------------------------------------------\n')
fprintf(['ClosedPut' fmtstring '\n'], EuroPut40(1,:))
fprintf(['MSPut    ' fmtstring '\n'], EuroPut40(2,:))
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('K = 45\n')
fprintf('-----------------------------------------------------------------------------------\n')
fprintf(['ClosedPut' fmtstring '\n'], EuroPut45(1,:))
fprintf(['MSPut    ' fmtstring '\n'], EuroPut45(2,:))
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('EuroPut = Black Scholes closed form European put price\n')
fprintf('MSPut   = Fifth-order Medvedev and Scaillet (2010) American put approximation\n')
fprintf('\n')
