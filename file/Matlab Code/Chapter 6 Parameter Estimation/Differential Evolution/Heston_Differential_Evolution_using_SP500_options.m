% Parameter estimation using loss functions and Fractional FFT
% Uses S&P500 options on April 13, 2012

clc; clear;

% Black Scholes call
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));

%% Load the SP500 data
MktIV = [ ...
   27.8000   26.3800   25.3200   25.1800
   24.7700   24.0200   23.6400   23.6900
   21.8600   21.5800   22.0300   22.3900
   18.7800   19.3000   20.4700   20.9800
   15.7200   17.1200   18.9400   19.7000
   13.3400   15.1700   17.4800   18.4900
   13.2300   13.7300   16.1800   17.3600]./100;
T = [45 98 261 348]./365;
K = [120   125   130   135   140   145   150];
[NK NT] = size(MktIV);
rf = 0.0010;
q  = 0.0068;
S = 137.14;
PutCall = repmat('C',NK,NT);


%% Use ITM Calls and Puts in the parameter estimation
for k=1:NK
    for t=1:NT
        MktPrice(k,t) = BSC(S,K(k),rf,q,MktIV(k,t),T(t));
    end
end


%% Parameter Estimation settings
% Initial values for kappa,theta,sigma,v0,rho
start = [9.0 0.05 0.3 0.05 -0.8];

% Specify the objective function and method to obtain the option price
ObjFun = 1;
Method = 1;

% Gauss Laguerre Weights and abscissas
[x w] = GenerateGaussLaguerre(32);
trap = 1;

% Settings for the Bisection method to find the Model IV
a = .001;
b = 3;
Tol = 1e-7;
MaxIter = 1000;

% Estimation bounds
e = 1e-5;
kappaL = e;       kappaU = 20;
thetaL = e;       thetaU = 2;
sigmaL = e;       sigmaU = 2;
v0L    = e;       v0U    = 2;
rhoL   = -.999;   rhoU   = .999;
Hi = [kappaU thetaU sigmaU v0U rhoU];
Lo = [kappaL thetaL sigmaL v0L rhoL];

% Gauss-Laguerre integration abscissas and weights
[x w] = GenerateGaussLaguerre(32);


%% Parameter Estimates using the FRFT
N = 2^7;
uplimit = 25;
eta = 0.25;
alpha = 1.75;
rule = 'T';
K1 = K(1)-1;
tic
FTparam = fmincon(@(p) HestonObjFunFRFT(p,S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule),start,[],[],[],[],Lo,Hi);
t1 = toc;


%% Differential Evolution algorithm
% Number of members in the population and number of generations
NP = 75;
NG = 200;

% Scale factor and crossover ratio
F  = 0.8;
CR = 0.5;

% Find the Differential Evolution parameters
tic
DEparam = HestonDE(NG,NP,CR,F,Hi,Lo,S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule);
t2 = toc;

%% Generate prices and parameters using both sets of parameters
[x w] = GenerateGaussLaguerre(32);
a = 0.001;
b = 5;
Tol = 1e-5;
MaxIter = 5000;
Sum1 = 0;
Sum2 = 0;

for k=1:NK
	for t=1:NT
		PriceFT(k,t)  = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,FTparam, trap,x,w);
		PriceDE(k,t)  = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,DEparam ,trap,x,w);
        IVFT(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,PriceFT(k,t),Tol,MaxIter);
        if IVFT(k,t) ~= -1;
            error1(k,t) = (IVFT(k,t) - MktIV(k,t))^2;
        else
            IVFT(k,t) = NaN;
            Sum1 = Sum1+1;
        end
		IVDE(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,PriceDE(k,t),Tol,MaxIter);
        if IVDE(k,t) ~= -1;
            error2(k,t) = (IVDE(k,t) - MktIV(k,t))^2;
        else
            IVDE(k,t) = NaN;
            Sum2 = Sum2+1;
        end
	end
end
ErrorFT  = sum(sum(error1)) / (NK*NT - Sum1);
ErrorDE =  sum(sum(error2)) / (NK*NT - Sum2);


%% Display the results
clc;
fmtstring = [repmat('%10.4f',1,6) '%12.2d'];
fprintf('--------------------------------------------------------------------------------\n')
fprintf('Method      kappa     theta     sigma      v0        rho      time      IVMSE\n')
fprintf('--------------------------------------------------------------------------------\n')
fprintf(['FRFT    ' fmtstring '\n'],FTparam,t1,ErrorFT)
fprintf(['DE      ' fmtstring '\n'],DEparam,t2,ErrorDE)
fprintf('--------------------------------------------------------------------------------\n')



%% Plot the results
for t=1:NT
	subplot(2,2,t)
	plot(K,MktIV(:,t),'kx:',K,IVFT(:,t),'kx-',K,IVDE(:,t),'ko-')
	legend('Market','Objective Function', 'Differential Evolution')
end

