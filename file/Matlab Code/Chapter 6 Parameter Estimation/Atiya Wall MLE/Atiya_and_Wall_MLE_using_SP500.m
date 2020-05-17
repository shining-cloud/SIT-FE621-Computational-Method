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
rf = 0.0033;
q  = 0.0063;
dt = 1/252;
Spot = xlsread('SPY Historical Prices.xls', 'table','G2');

% Generate market prices from the quoted implied vols
for t=1:NT
    for k=1:NK
        MktPrice(k,t) = BSC(Spot,K(k),rf,q,MktIV(k,t),T(t));
    end
end

%% Parameter estimates by minimizing the negative Likelihood.m function
% Select method : 1 = Likelihood, 2 = Log-Likelihood.  Set the options.
Lmethod = 2;
options = optimset('MaxFunEvals',1e5,'MaxIter',1e5);

% Bounds on the parameter estimates and starting values
e = 1e-5;
% kappa theta sigma v0 rho
lb = [e   e  e  e  -0.999];  % Lower bound on the estimates
ub = [30  3  3  2   0.999];  % Upper bound on the estimates
start = [5 0.1 0.5 0.1 -0.8];

% Price calculation method and objective function
Method = 1;
ObjFun = 4;
trap = 1;

% Settings for the bisection algorithm
a = 0.01;
b = 3.0;
Tol = 1e-5;
MaxIter = 1000;

% Gauss Laguerre abscissas and weights
[X W] = GenerateGaussLaguerre(32);


%% Estimate parameters using sets of stock prices for each maturity
for t=1:NT
    % Historical stock prices.  Use D days before May 7, 2010.
    D(t) = T(t)*365;
    S = xlsread('SPY Historical Prices.xls', 'table',['G2:G' num2str(D(t)+1)]);

    % Put oldest prices first and calculate log prices
    S = flipud(S);
    x = log(S);

    % Obtain the parameter estimates using option prices and loss functions
    [true(t,:) feval(t)] = fmincon(@(p) HestonObjFun(p,Spot,rf,q,MktPrice(:,t),K,T(t),PutCall,MktIV(:,t),X,W,trap,ObjFun,Method,a,b,Tol,MaxIter),start,[],[],[],[],lb,ub,[],options);

    % Obtain the parameters estimates using Atiya and Wall MLE
    param(t,:) = fmincon(@(p) LikelihoodAW(p,x,rf,q,dt,Lmethod),start,[],[],[],[],lb,ub,[],options);
end

%% Generate the implied volatilities from the Atiya and Wall parameters
for t=1:NT
    N(t) = 0;
    for k=1:NK
        PriceAW(k,t) = HestonPriceGaussLaguerre('C',Spot,K(k),T(t),rf,q,param(t,:),trap,X,W);
        IVAW(k,t) = BisecBSIV('C',Spot,K(k),rf,q,T(t),a,b,PriceAW(k,t),Tol,MaxIter);
        if IVAW(k,t) ~= -1;
            N(t) = N(t)+1;
            Error(k,t) = abs(IVAW(k,t) - MktIV(k,t))/MktIV(k,t) * 100;
        else
            IVAW(k,t) = NaN;
            Error(k,t) = 0;
        end
    end
    IVRMSE(t) = sum(Error(:,t))/N(t);
end

%% Display the true parameter values and compare to the estimates
clc;
fprintf('Estimates                       kappa     theta     sigma     v0        rho     IVRMSE\n')
fprintf('----------------------------------------------------------------------------------------\n')
fmtstring = repmat('%10.4f',1,5);
for t=1:NT
    if Lmethod==1
        fprintf('Log-Likelihood %5.3d for estimates using %3.0d days\n',log(-feval(t)),D(t));
    elseif Lmethod==2
        fprintf('Log-Likelihood %5.3d for estimates using %3.0d days\n',-feval(t),D(t));
    end
    fprintf(['Atiya-Wall MLE Estimates    ' fmtstring ' %10.4f\n'], param(t,:), IVRMSE(t))
    fprintf(['Loss Function Estimates     ' fmtstring '\n'], true(t,:))
fprintf('----------------------------------------------------------------------------------------\n')
end

%% Plot the implied volatilities
for t=1:NT
	subplot(2,2,t)
	plot(K,MktIV(:,t),'bx-',K,IVAW(:,t),'ro-')
	title(['Maturity ' num2str(T(t)*365) ' days'])
	legend('Market IV', 'Heston IV')
end

