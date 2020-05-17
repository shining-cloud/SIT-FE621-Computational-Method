% Maximum Likelihood Estimates of the Heston parameters,
% using the method of A. Atiya and S. Wall

clc; clear;

% Function for the Black Scholes call
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));

%% Read in the data 
% Implied volatilities from calls
MktIV = [ ...
   27.8000   26.3800   25.3200   25.1800
   24.7700   24.0200   23.6400   23.6900
   21.8600   21.5800   22.0300   22.3900
   18.7800   19.3000   20.4700   20.9800
   15.7200   17.1200   18.9400   19.7000
   13.3400   15.1700   17.4800   18.4900
   13.2300   13.7300   16.1800   17.3600]./100;
[NK NT] = size(MktIV);
PutCall = repmat('C',NK,NT);

% Maturities, strikes, rates, time increment, Spot
T = [45 98 261 348]./365;
K = [120   125   130   135   140   145   150];
rf = 0.0010;
q  = 0.0068;
dt = 1/252;
Spot = xlsread('SPY Historical Prices.xls', 'table','G2');

% Generate market prices from the quoted implied vols
for t=1:NT
    for k=1:NK
        MktPrice(k,t) = BSC(Spot,K(k),rf,q,MktIV(k,t),T(t));
    end
end

%% Loss function parameter estimates

% Optimization options
options = optimset('MaxFunEvals',1e5,'MaxIter',1e5);

% kappa theta sigma v0 rho
e = 1e-5;
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [20  2  2  3  .999];  % Upper bound on the estimates
start = [9.0 0.05 0.3 0.05 -0.8];

% Price calculation method and objective function
Method = 1;
ObjFun = 1;
trap = 1;

% Settings for the bisection algorithm
a = 0.01;
b = 3.0;
Tol = 1e-5;
MaxIter = 1000;

% Gauss Laguerre abscissas and weights
[X W] = GenerateGaussLaguerre(32);

% Obtain the parameter estimates using option prices and loss functions
[true feval] = fmincon(@(p) HestonObjFun(p,Spot,rf,q,MktPrice,K,T,PutCall,MktIV,X,W,trap,ObjFun,Method,a,b,Tol,MaxIter),start,[],[],[],[],lb,ub,[],options);

%% MLE parameters using sets of stock prices for each maturity

% Select method : 1 = Likelihood, 2 = Log-Likelihood.  Set the options.
Lmethod = 2;

% Choose the number of days to use in the MLE. Run code again with t=2,t=3,t=4.
t = 2;

% Historical stock prices.  Use D days before May 7, 2010.
D = T(t)*365;
S = xlsread('SPY Historical Prices.xls', 'table',['G2:G' num2str(D+1)]);

% Put oldest prices first and calculate log prices
S = flipud(S);
x = log(S);

% Obtain the parameters estimates using Atiya and Wall MLE
param = fmincon(@(p) LikelihoodAW(p,x,rf,q,dt,Lmethod),start,[],[],[],[],lb,ub,[],options);

%% Generate the implied volatilities from the Atiya and Wall parameters
for t=1:NT
    N(t) = 0;
    for k=1:NK
        PriceAW(k,t) = HestonPriceGaussLaguerre('C',Spot,K(k),T(t),rf,q,param,trap,X,W);
        IVAW(k,t) = BisecBSIV('C',Spot,K(k),rf,q,T(t),a,b,PriceAW(k,t),Tol,MaxIter);
        if IVAW(k,t) ~= -1;
            N(t) = N(t)+1;
            Error(k,t) = (IVAW(k,t) - MktIV(k,t))^2;
        else
            IVAW(k,t) = NaN;
            Error(k,t) = 0;
        end
    end
    IVRMSE(t) = sum(Error(:,t))/N(t);
end
IVRMSE = sum(IVRMSE)/NK;

%% Display the true parameter values and compare to the estimates
clc;
fprintf('Estimates                       kappa     theta     sigma     v0        rho     IVRMSE\n')
fprintf('----------------------------------------------------------------------------------------\n')
fmtstring = repmat('%10.4f',1,5);
if Lmethod==1
    fprintf('Log-Likelihood %5.3d for estimates using %3.0d days\n',log(-feval),D);
elseif Lmethod==2
    fprintf('Log-Likelihood %5.3d for estimates using %3.0d days\n',-feval,D);
end
fprintf(['Atiya-Wall MLE Estimates    ' fmtstring ' %10.4f\n'], param, IVRMSE)
fprintf(['Loss Function Estimates     ' fmtstring '\n'], true)
fprintf('----------------------------------------------------------------------------------------\n')

%% Plot the implied volatilities
for t=1:NT
	subplot(2,2,t)
	plot(K,MktIV(:,t),'bx-',K,IVAW(:,t),'ro-')
	title(['Maturity ' num2str(T(t)*365) ' days'])
	legend('Market IV', 'AW-Heston IV')
end

