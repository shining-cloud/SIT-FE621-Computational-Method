clc; clear;

K = [1100:5:1225];
S  = 1164.97;
rf = 0.0;
q  = 0.0071;
T = 15/365;
trap = 1;

NK = length(K);
kappa =  3.0002;
theta =  0.0633;
sigma =  0.5015;
v0    =  0.1044;
rho   = -0.9008;
param = [kappa theta sigma v0 rho];

a = 0.01;
b = 2.0;
Tol = 1e-5;
MaxIter = 1000;

%% Estimate of fair vol by Heston
Hvar = (v0-theta)*(1-exp(-kappa*T))/kappa/T + theta;
Hvol = sqrt(Hvar).*100;


%% Demeterfi et al variance swap fair strike 
% Generate calls and puts 
[x w] = GenerateGaussLaguerre(32);

for k=1:NK
    Call(k) = HestonPriceGaussLaguerre('C',S,K(k),T,rf,q,param,trap,x,w);
    Put(k)  = HestonPriceGaussLaguerre('P',S,K(k),T,rf,q,param,trap,x,w);
    CallIV(k) = BisecBSIV('C',S,K(k),rf,q,T,a,b,Call(k),Tol,MaxIter);
    PutIV(k)  = BisecBSIV('P',S,K(k),rf,q,T,a,b,Put(k),Tol,MaxIter);
end;


%% Select the IV for OTM calls and puts and the number of interpolation points
ATM = find(K==round(S));
NI = 2000;

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

%% Calculate variance swap fair strike
Kvar = VarianceSwap(KCI,CallVI,KPI,PutVI,S,T,rf,q);
Kvol = sqrt(Kvar)*100;


%% Display the results
disp(['Using ' num2str(NI) ' interpolation points'])
disp(['Estimate of fair volatility using replication       ' num2str(Kvar,'%10.8f')])
disp(['Estimate of fair volatility from Heston parameters  ' num2str(Hvar,'%10.8f')])


