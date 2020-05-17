% Parameter estimation using loss functions and Fractional FFT
% Uses IBM options on May 7, 2010

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

%% Use Calls in the parameter estimation
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
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [20  2  2  2  .999];  % Upper bound on the estimates

% Gauss-Laguerre integration abscissas and weights
[x w] = GenerateGaussLaguerre(32);

%% Parameter Estimates using the regular objective function
tic
param = fmincon(@(p) HestonObjFun(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,Method,a,b,Tol,MaxIter),start,[],[],[],[],lb,ub);
t1 = toc;

%% Parameter Estimates using the SVC objective function, Heston CF
CF = 'Heston';
tic
paramSVC1 = fmincon(@(p) HestonObjFunSVC(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,CF),start,[],[],[],[],lb,ub);
t2 = toc;

%% Parameter Estimates using the SVC objective function, Attari CF
CF = 'Attari';
tic
paramSVC2 = fmincon(@(p) HestonObjFunSVC(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,CF),start,[],[],[],[],lb,ub);
t3 = toc;

%% Parameter Estimates using the FRFT
N = 2^7;
uplimit = 25;
eta = 0.25;
alpha = 1.75;
rule = 'T';
K1 = K(1)-1;
tic
paramFRFT = fmincon(@(p) HestonObjFunFRFT(p,S,K1,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,a,b,Tol,MaxIter,trap,N,uplimit,eta,alpha,rule),start,[],[],[],[],lb,ub);
t4 = toc;


%% Fit the model implied volatilities from the model prices
% "Little Trap" formulation for the Heston characteristic function
trap = 1;
Sum0 = 0;
Sum1 = 0;
Sum2 = 0;
SumF = 0;
for k=1:NK
	for t=1:NT
 		OrdPrice(k,t)  = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,param,    trap,x,w);
 		SVCPrice1(k,t) = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,paramSVC1,trap,x,w);
 		SVCPrice2(k,t) = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,paramSVC2,trap,x,w);
 		FRFTPrice(k,t) = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,paramFRFT,trap,x,w);
        OrdIV(k,t)  = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,OrdPrice(k,t), Tol,MaxIter);
        SVCIV1(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,SVCPrice1(k,t),Tol,MaxIter);
		SVCIV2(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,SVCPrice2(k,t),Tol,MaxIter);
        FRFTIV(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,FRFTPrice(k,t),Tol,MaxIter);
        if OrdIV(k,t) ~= -1;
            error0(k,t) = (OrdIV(k,t) - MktIV(k,t))^2;
        else
            OrdIV(k,t) = NaN;
            Sum0 = Sum0+1;
        end
        if SVCIV1(k,t) ~= -1;
            error1(k,t) = (SVCIV1(k,t) - MktIV(k,t))^2;
        else
            SVCIV1(k,t) = NaN;
            Sum1 = Sum1+1;
        end
        if SVCIV2(k,t) ~= -1;
            error2(k,t) = (SVCIV2(k,t) - MktIV(k,t))^2;
        else
            SVCIV2(k,t) = NaN;
            Sum2 = Sum2+1;
        end
        if FRFTIV(k,t) ~= -1;
            errorF(k,t) = (FRFTIV(k,t) - MktIV(k,t))^2;
        else
            FRFTIV(k,t) = NaN;
            SumF = SumF+1;
        end
	end
end
ErrorOrd   = sum(sum(error0)) / (NK*NT - Sum0);
ErrorSVC1  = sum(sum(error1)) / (NK*NT - Sum1);
ErrorSVC2  = sum(sum(error2)) / (NK*NT - Sum2);
ErrorFRFT  = sum(sum(errorF)) / (NK*NT - SumF);

 
%% Plot the implied volatilities
for t=1:NT
	subplot(2,2,t)
	plot(K,MktIV(:,t),'k:',K,OrdIV(:,t), K,SVCIV1(:,t),'kx-',K,SVCIV2(:,t),'kv-',K,FRFTIV(:,t),'ko-')
	title(['Maturity ' num2str(T(t)*365.25) ' days'])
	legend('Market Implied Volatility', 'Ordinary IV Heston', 'SVC IV Heston', 'SVC IV Attari', 'FRFT IV')
	xlabel('Strike Price')
	ylabel('Implied Volatility')
end

%% Display the results
clc
fmtstring = [repmat('%10.4f',1,6) '%12.2d'];
fprintf('----------------------------------------------------------------------------------------------\n')
fprintf('Objective Function       kappa     theta     sigma     v0        rho      time      IVMSE\n')
fprintf('----------------------------------------------------------------------------------------------\n')
fprintf(['Ordinary             ' fmtstring '\n'], param,    t1,ErrorOrd)
fprintf(['SVC-Heston           ' fmtstring '\n'], paramSVC1,t2,ErrorSVC1)
fprintf(['SVC-Attari           ' fmtstring '\n'], paramSVC2,t3,ErrorSVC2)
fprintf(['Fractional FFT       ' fmtstring '\n'], paramFRFT,t4,ErrorFRFT)
fprintf('----------------------------------------------------------------------------------------------\n')

