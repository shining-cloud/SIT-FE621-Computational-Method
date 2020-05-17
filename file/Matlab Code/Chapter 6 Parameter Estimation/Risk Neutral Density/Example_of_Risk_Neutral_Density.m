clc; clear;

% Maturities 
T = [45 98 261 348]./365;
NT = length(T);

% Option settings
S0 = 137.14;
rf = 0.0010;
q  = 0.0068;

%% Heston parameters from SP data, Table 6.1
% MSE parameters
% kappa =  1.9214;
% theta =  0.0904;
% sigma =  1.0193;
% v0    =  0.0344;
% rho   = -0.7799;
% lambda = 0;

% RMSE parameters
kappa =  8.9931;
theta =  0.0571;
sigma =  2.0000;
v0    =  0.0405;
rho   = -0.7899;
lambda = 0;

% IVRMSE parameters
% kappa =  9.0000;
% theta =  0.0420;
% sigma =  0.3044;
% v0    =  0.0405;
% rho   = -0.8038;
% lambda = 0;

% CHJ (2009) parameters
% kappa =  8.9655;
% theta =  0.0602;
% sigma =  1.7142;
% v0    =  0.7913;
% rho   = -0.9921;
% lambda = 0;

param = [kappa theta sigma v0 rho];
trap = 1;

% Gauss-Legendre abscissas and weights
[xGLe wGLe] = GenerateGaussLegendre(32);

% Multi-domain integration rule settings
lo = 1e-10;
hi = 500;
Ndomains  = 50 ;
dA = (hi - lo)/Ndomains;
A = [lo:dA:hi];
tol = 1e-8;

%% Initialize the cells
Call   = cell(NT,1);
Strike = cell(NT,1);
Domain = cell(NT,1);
Npoint = cell(NT,1);
RND = cell(NT,1);
K   = cell(NT,1);

%% Construct the Risk Neutral Densities

% Strike increment and ranges
dK = 0.5;
% Use these strikes to plot the RND
K{1} = 90:dK:170;
K{2} = 80:dK:180;
K{3} = 40:dK:220;
K{4} = 20:dK:240;

% Use these strikes to recover the market call pri
% K{1} = 70:dK:190;
% K{2} = 60:dK:200;
% K{3} = 20:dK:240;
% K{4} = 10:dK:260;

% Extract the RND, integration domain, area, negative values
for t=1:NT;
    NK = length(K{t});
    for k=1:NK;
        [Call{t}(k) L H(k) N{t}(k)] = HestonPriceGaussLegendreMD('C',S0,K{t}(k),T(t),rf,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLe,wGLe,A,tol);
    end
    Domain{t} = max(H);
    [RND{t} ST{t}] = ExtractRND(K{t},Call{t});
    Area(t) = trapz(RND{t})*dK;
    Zero(t) = length(find(RND{t}<0));
end

%% Output the results
fprintf('Risk Neutral Density estimation \n');
fprintf('Maturity       Area   NegValues  IntLimit  IntPoints\n');
fprintf('----------------------------------------------------\n');
for t=1:NT
    fprintf(' %5.0f    %10.4f  %6.0f %10.0f %8.0f \n',T(t)*365,Area(t),Zero(t),Domain{t},max(N{t}));
end
fprintf('----------------------------------------------------\n');

 
%% Plot the RNDs
plot(ST{4},RND{4},'c-',ST{3},RND{3},'b-',ST{2},RND{2},'r-',ST{1},RND{1},'k')
legend([num2str(T(4)*365) '-day maturity'],...
       [num2str(T(3)*365) '-day maturity'],...
       [num2str(T(2)*365) '-day maturity'],...
       [num2str(T(1)*365) '-day maturity'])
xlabel('Terminal Stock Price');
ylabel('RND');
xlim([60 200]);

%% Recover the calls
% Market prices and strikes from data used to generate Table 6.1
MktPrice = [...
   17.5399   18.4781   21.1350   22.3635;
   12.8889   14.1227   17.1876   18.5477;
    8.5359   10.0404   13.5115   15.0476;
    4.6903    6.4381   10.1650   11.7645;
    1.7960    3.5134    7.2126    8.8694;
    0.3665    1.5057    4.7473    6.3808;
    0.0654    0.4821    2.8619    4.3398];
MktStrike = [120   125   130   135   140   145   150];
NK = length(MktStrike);

for t=1:NT
    for k=1:NK
        Payoff = max(ST{t} - MktStrike(k), 0);
        RNDCall(k,t) = trapz(Payoff.*RND{t})*dK*exp(-rf*T(t));
    end
end
error = sum(sum(MktPrice - RNDCall));

%% Print the results
fprintf('Market prices from data used to generate Table 6.1\n')
fprintf('Strike   45-day     98-day     261-day    348-day \n')
fprintf('--------------------------------------------------\n')
for k=1:7
    fprintf('%5.0f %10.4f %10.4f %10.4f %10.4f \n',MktStrike(k),MktPrice(k,1),MktPrice(k,2),MktPrice(k,3),MktPrice(k,4));
end
fprintf('--------------------------------------------------\n')

fprintf(' \n');
fprintf('Recovered call prices  \n')
fprintf('Strike   45-day     98-day     261-day    348-day \n')
fprintf('--------------------------------------------------\n')
for k=1:7
    fprintf('%5.0f %10.4f %10.4f %10.4f %10.4f \n',MktStrike(k),RNDCall(k,1),RNDCall(k,2),RNDCall(k,3),RNDCall(k,4));
end
fprintf('--------------------------------------------------\n')
fprintf('Sum of absolute errors %5.4f \n',error)
fprintf('--------------------------------------------------\n')

