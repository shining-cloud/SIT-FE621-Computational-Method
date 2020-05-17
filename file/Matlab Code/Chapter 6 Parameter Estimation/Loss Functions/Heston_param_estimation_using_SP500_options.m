% Parameter estimation using loss functions.
% Uses SP500 options on April 14, 2012

clc; clear;

% Black Scholes call and put
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));

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


%% Use ITM Calls in the parameter estimation
for k=1:NK
    for t=1:NT
        MktPrice(k,t) = BSC(S,K(k),rf,q,MktIV(k,t),T(t));
    end
end

%% Parameter Estimation
% Initial values for kappa,theta,sigma,v0,rho
start = [9.0 0.05 0.3 0.05 -0.8];

% Specify the objective function
% 1 = MSE | 2 = RMSE | 3 = IVMSE | 4 = Christoffersen, Heston, Jacobs
ObjFun = 1;

% Specify the method to obtain option prices
%   1 = Heston C.F. with Gauss Laguerre integration
%   2 = Lewis vol-of-vol Series II expansion
%   3 = Lewis vol-of-vol Series I expansion
Method = 1;

% Gauss Laguerre Weights and abscissas
[x w] = GenerateGaussLaguerre(32);
trap = 1;

%% Find the parameter estimates

% Settings for the Bisection method to find the Model IV
a = .001;
b = 3;
Tol = 1e-7;
MaxIter = 1000;

e = 1e-5;
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [20  2  2  2  .999];  % Upper bound on the estimates
tic
[param feval] = fmincon(@(p) HestonObjFun(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,Method,a,b,Tol,MaxIter),start,[],[],[],[],lb,ub)
t1 = toc;

% Parameter estimates
kappa  = param(1);
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);
lambda = 0;

%% Fit the model implied volatilities from the model prices
% "Little Trap" formulation for the Heston characteristic function
trap = 1;

% Integration rule settings
method = 3;
a = 1e-10;
b = 250;
N = 15000;

error = zeros(NK,NT);
Sum = 0;
for k=1:NK
	for t=1:NT
 		ModelPrice(k,t) = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,param,trap,x,w);
		ModelIV(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
        if ModelIV(k,t) ~= -1
            error(k,t) = (ModelIV(k,t) - MktIV(k,t))^2;
        else
            Sum = Sum + 1;
        end
	end
end
ErrorIV = sum(sum(error)) / (NK*NT - Sum);

%% Display the results
clc;
fprintf('Loss function estimation using loss function %1.0f \n',ObjFun);
fprintf('-----------------------------------------------------------------------\n');
fprintf('  kappa      theta      sigma       v0        rho    IVMSE     EstTime\n');
fprintf('-----------------------------------------------------------------------\n');
fprintf('%7.4f %10.4f %10.4f %10.4f %10.4f %10.2e %8.3f\n',param,ErrorIV,t1);

%% Plot the implied volatilities
for t=1:NT
	subplot(2,2,t)
	plot(K,MktIV(:,t),'kx-',K,ModelIV(:,t),'ro-')
	title(['Maturity ' num2str(T(t)*365) ' days'])
	legend('Market IV', 'Heston IV')
end

%% Plot the total volatilities
% subplot(1,1,1)
% plot(K,ModelIV(:,1).*sqrt(T(1)),'ko-',...
%      K,ModelIV(:,2).*sqrt(T(2)),'ro-',...
%      K,ModelIV(:,3).*sqrt(T(3)),'bo-',...
%      K,ModelIV(:,4).*sqrt(T(4)),'co-')

