clc; clear;

%% Read in the raw data.  Use Calls or Puts
Call = xlsread('SPX_options.xls', 'Table', 'B67:E92');   % Calls
Put  = xlsread('SPX_options.xls', 'Table', 'M67:P92');   % Puts

% Read in the maturities and strikes
K = xlsread('SPX_options.xls', 'Table', 'A67:A92');
T = xlsread('SPX_options.xls', 'Table', 'M1:P1');

%% Select a single maturity if desired
t = 2;
Call = Call(:,t);
Put  = Put(:,t);
T = T(t);
[NK,NT] = size(Call);

%% Shimko's regression to get implied dividend yield and risk-free rate.
% intercept is an estimate of S*D(t,T) and
% slope is an estimate of B(t,T).
% Spot price
S  = 1164.97;

for i=1:NT
	Y = Call(:,i) - Put(:,i);
	B = regstats(Y,K,'linear');
	int(i)   = B.beta(1);
	slope(i) = B.beta(2);
	riskfree(i) = -log(-slope(i))/T(i);
	riskfree(i) = max(riskfree(i),0);
	dividend(i) = log(S/int(i))/T(i);
	dividend(i) = max(dividend(i),0);
end

% Risk free rate and dividend yield
rf = mean(riskfree);
q  = mean(dividend);

% Implied volatilities
a = 0.001;
b = 10;
Tol = 1e-10;
MaxIter = 1000;

for k=1:NK
    CallIV(k) = BisecBSIV('C',S,K(k),rf,q,T,a,b,Call(k),Tol,MaxIter);
     PutIV(k) = BisecBSIV('P',S,K(k),rf,q,T,a,b,Put(k),Tol,MaxIter);
end


%% Replication algorithm on market quotes

% Find the IV for OTM calls and puts and select # of interpolation points
ATM = find(K==round(S));
NI = 1000;

%% Select OTM strikes for puts, create fine grid, and interpolate the IV
KC = K(ATM:NK);
CallV = CallIV(ATM:NK);
dKC = (KC(end) - KC(1))/NI;
KCI = [KC(1):dKC:KC(end)];
CallVI = interp1(KC,CallV,KCI,'linear');
 
%% Select OTM strikes for puts, create fine grid, and interpolate the IV
KP = K(1:ATM);
PutV = PutIV(1:ATM);
dKP = (KP(end) - KP(1))/NI;
KPI = [KP(1):dKP:KP(end)];
PutVI = interp1(KP,PutV,KPI,'linear');


%% Find the variance swap strike
Kvar = VarianceSwap(KCI,CallVI,KPI,PutVI,S,T,rf,q);


%% Use calls or puts to estimate the Heston parameters
Choice = 'C';

if strcmp(Choice,'C')
    PutCall = repmat('C',NK,1);
    MktPrice = Call;
    MktIV = CallIV';
else
    PutCall = repmat('P',NK,1);
    MktPrice = Put;
    MktIV = PutIV';
end

%% Parameter Estimation
% kappa theta sigma v0 rho
start = [3 .06 .5 .05 -0.9];

% Specify the objective function
ObjFun = 4;

% Weights and abscissas
[x w] = GenerateGaussLaguerre(32);

% Find the parameter estimates
e = 1e-3;
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [10 3  3  3 .999];  % Upper bound on the estimates
CF = 'Heston';
trap = 1;
[param fe] = fmincon(@(beta) HestonObjFunSVC(beta,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,CF),start,[],[],[],[],lb,ub);

% Parameter estimates
kappa  = param(1); 
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);
lambda = 0;

%% Heston variance swap
Hvar = (v0-theta)*(1-exp(-kappa*T))/kappa/T + theta;


%% Display the results
fprintf('Heston parameter estimates\n')
fprintf('kappa %5.4f\n',kappa)
fprintf('theta %5.4f\n',theta)
fprintf('sigma %5.4f\n',sigma)
fprintf('v0    %5.4f\n',v0)
fprintf('rho   %5.4f\n',rho)
fprintf('------------------------------------------------------------------\n')
fprintf('Estimate of fair volatility using replication         %10.5f\n', Kvar)
fprintf('Estimate of fair volatility using Heston parameters   %10.5f\n', Hvar)
fprintf('------------------------------------------------------------------\n')

