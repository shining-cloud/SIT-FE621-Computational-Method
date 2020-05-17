% Parameter estimation starting values 
% using method of Gauthier and Rivaille

clc; clear;

%% Load the data: 4 maturities, 7 strikes
Calls = [...
   18.2500   20.5650   21.6200   24.4600
   13.9150   16.6600   17.8550   21.0550
    9.8600   13.0350   14.4200   17.9450
    6.2950    9.7500   11.2100   15.0000
    3.4150    6.8700    8.3900   12.2900
    1.4500    4.4800    5.9850    9.8750
    0.4600    2.6700    4.0300    7.7500];
Puts = [...
    1.6200    4.9650    6.6200   11.3000
    2.3250    6.1400    7.9700   12.9800
    3.3450    7.6000    9.5500   14.8900
    4.8400    9.3800   11.3950   16.9850
    7.0100   11.6150   13.6250   19.1850
   10.0550   14.3050   16.3200   21.5500
   14.0850   17.5550   19.3800   24.4350];
PutIV = [ ...
   26.3500   25.8200   25.9400   25.9500
   24.1900   24.3000   24.6300   24.9700
   22.0200   22.8300   23.3200   24.0600
   19.9200   21.3900   21.9900   23.1300
   17.9300   20.1200   20.8000   22.0400
   16.2100   18.9200   19.8400   20.8500
   15.5600   17.9800   18.8700   20.0500]./100;
T = [0.2685   0.7151    0.9534    1.7644];
Strike = [120   125   130   135   140   145   150];
[NK NT] = size(Calls);
S = 137.14;


%% Settings for the Gauthier method to find starting values
% Select the maturity on which to implement the method
t = 2;

% Select arbitrary starting values for kappa, theta, and v0
kappa0 = 10;
theta0 = 0.1;
v00    = 0.01;
sigma0 = 0.3;
rho0   = -0.7;

% Select the two puts and strikes on which to implement the Gauthier method
I = [3 7];
K = Strike(I);
tau = T(t);
Put1 = Puts(I(1),t);
Put2 = Puts(I(2),t);

% Shimko's regression to get implied dividend yield and risk-free rate.
Y = Calls(:,t) - Puts(:,t);
B = regstats(Y,Strike,'linear');
int   = B.beta(1);
slope = B.beta(2);
rf = -log(-slope)/T(t);
rf = max(rf,0);
q  = log(S/int)/T(t);
q = max(q,0);

%% Get the Gauthier starting values
method = 2;
[sigma rho] = GetGauthierValues(kappa0,theta0,v00,S,K(1),K(2),Put1,Put2,tau,rf,q,method)


%% Use the starting values in the optimization
[x w] = GenerateGaussLaguerre(32);
trap = 1;

% Use Puts in the parameter estimation
MktIV    = PutIV(:,t);
MktPrice = Puts(:,t);
[NK,NT] = size(MktPrice);
PutCall = repmat('P',NK,NT);

% Settings for the Bisection method to find the Model IV
a = .01;
b = 2;
Tol = 1e-5;
MaxIter = 1000;
e = 1e-5;
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [100 10 10 10 .999];  % Upper bound on the estimates
ObjFun = 1;  % Select the loss function
Method = 1;  % Select the method to obtain Heston prices

% Parameter Estimates using arbitrary starting values
start1 = [kappa0 theta0 sigma0 v00 rho0];
tic
param1 = fmincon(@(p) HestonObjFunSVC(p,S,rf,q,MktPrice,Strike,tau,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,'Heston'),start1,[],[],[],[],lb,ub);
t1 = toc;

% Parameter Estimates using Gauthier starting values
start2 = [kappa0 theta0 sigma v00 rho];
tic
param2 = fmincon(@(p) HestonObjFunSVC(p,S,rf,q,MktPrice,Strike,tau,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,'Heston'),start2,[],[],[],[],lb,ub);
t2 = toc;

%% Find the implied volatilities for each set of parameter estimates
for k=1:length(Strike)
    Price1(k) = HestonPriceGaussLaguerre(PutCall(k),S,Strike(k),tau,rf,q,param1,trap,x,w);
    Price2(k) = HestonPriceGaussLaguerre(PutCall(k),S,Strike(k),tau,rf,q,param2,trap,x,w);
    IV1(k) = BisecBSIV(PutCall(k),S,Strike(k),rf,q,tau,a,b,Price1(k),Tol,MaxIter);
    IV2(k) = BisecBSIV(PutCall(k),S,Strike(k),rf,q,tau,a,b,Price2(k),Tol,MaxIter);
    error1(k) = (IV1(k) - MktIV(k))^2;
    error2(k) = (IV2(k) - MktIV(k))^2;
end
IVMSE1 = sum(error1)/NK;
IVMSE2 = sum(error2)/NK;


%% Display the results
clc
fprintf('Gauthier Starting values: %7.4f %4.4f \n',sigma,rho)
fprintf(' \n');
fprintf('Starting Values     kappa    theta     sigma     v0      rho     time     IVMSE\n');
fprintf('---------------------------------------------------------------------------------\n');
fprintf('Arbitrary  %15.4f %8.4f %8.4f %8.4f %8.4f %8.4f %10.2e\n',param1(1),param1(2),param1(3),param1(4),param1(5),t1,IVMSE1);
fprintf('Gauthier   %15.4f %8.4f %8.4f %8.4f %8.4f %8.4f %10.2e\n',param2(1),param2(2),param2(3),param2(4),param2(5),t2,IVMSE2);

