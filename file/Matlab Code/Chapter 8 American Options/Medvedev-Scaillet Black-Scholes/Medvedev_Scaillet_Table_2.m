% Reproduces Table 2 in Medvedev and Scaillet (2010) for American puts
% under the Black Scholes model

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

start = 2;
options = optimset('LargeScale','off');


%% Find the tree and M-S prices for for Table 2
K = [35 40 45];
for k=1:3
    for i=1:9
        BSPut(k,i) = BSP(S,K(k),r,q,sigma(i),T(i));
        [MSPut(k,i) y(k,i) theta(k,i)] = MSPriceBS(S,K(k),T(i),sigma(i),r,q);
        BTPut(k,i) = TrinomialTree(S,K(k),r,q,sigma(i),T(i),N,PutCall,EuroAmer);
    end
end


%% Output the results
fmtstring = repmat('%8.5f',1,9);
fprintf('                           Table 2 of Medvedev and Scaillet (2010)\n');
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('                sigma = 0.2             sigma = 0.3             sigma = 0.4\n')
fprintf('           ----------------------  ----------------------  ------------------------\n');
fprintf('           T=1/12  T=1/3   T=7/12  T=1/12  T=1/3   T=7/12  T=1/12  T=1/3   T=7/12  \n');
fprintf('-----------------------------------------------------------------------------------\n')
for k=1:3
    fprintf(['K = ' num2str(K(k)) '\n'])
    fprintf(['EuroPut   ' fmtstring '\n'], BSPut(k,:))
    fprintf(['MSPut     ' fmtstring '\n'], MSPut(k,:))
    fprintf(['TreePut   ' fmtstring '\n'], BTPut(k,:))
    fprintf('-----------------------------------------------------------------------------------\n')
end
fprintf('EuroPut = Black Scholes closed form European put price\n')
fprintf('MSPut   = Fifth-order Medvedev and Scaillet (2010) American put approximation\n')
fprintf('TreePut = Trinomial tree American put with %3.0f steps\n',N);
fprintf('\n')

%% Output the barriers and moneyness
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('                sigma = 0.2             sigma = 0.3             sigma = 0.4\n')
fprintf('           ----------------------  ----------------------  ------------------------\n');
fprintf('           T=1/12  T=1/3   T=7/12  T=1/12  T=1/3   T=7/12  T=1/12  T=1/3   T=7/12  \n');
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('Barrier levels y \n');
fprintf(['K = 35   ' fmtstring '\n'], y(1,:))
fprintf(['K = 40   ' fmtstring '\n'], y(2,:))
fprintf(['K = 45   ' fmtstring '\n'], y(3,:))
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('Moneyness (theta) \n');
fprintf(['K = 35   ' fmtstring '\n'], theta(1,:))
fprintf(['K = 40   ' fmtstring '\n'], theta(2,:))
fprintf(['K = 45   ' fmtstring '\n'], theta(3,:))
fprintf('\n')
fprintf('-----------------------------------------------------------------------------------\n')

%% Check barriers are >= theta
y >= theta
