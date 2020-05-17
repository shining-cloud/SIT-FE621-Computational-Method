% Parameter estimation using loss functions.
% Uses DJIA options

clc; clear;

% Black Scholes call and put
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));

%% Load the DJIA data
MktIV = [...
    19.62	19.47	20.19	21.15;
    19.10	19.05	19.80	20.82;
    18.60	18.61	19.43	20.57;
    18.10	18.12	19.07	20.21;
    17.61	17.74	18.71	20.00;
    17.18	17.43	18.42	19.74;
    16.71	17.06	18.13	19.50;
    16.44	16.71	17.83	19.27;
    16.45	16.41	17.60	18.99;
    16.61	16.25	17.43	18.84;
    17.01	16.02	17.26	18.62;
    17.55	16.10	17.16	18.46;
    17.96	16.57	17.24	18.42]./100;
K = (124:136);
T = [37 72 135 226]./365;
[NK NT] = size(MktIV);
S = 129.14;
rf = 0;
q  = 0;
PutCall = repmat('P',NK,NT);


%% Use ITM Calls in the parameter estimation
for k=1:NK
    for t=1:NT
        MktPrice(k,t) = BSP(S,K(k),rf,q,MktIV(k,t),T(t));
    end
end

%% Parameter Estimation
% Initial values for kappa,theta,sigma,v0,rho
start = [9.0 0.05 0.3 0.05 -0.8];

% Specify the objective function
% 1 = MSE | 2 = RMSE | 3 = IVMSE | 4 = Christoffersen, Heston, Jacobs
ObjFun = 2;

% Specify the method to obtain option prices
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
e = 1e-2;
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [20  5  5  5  .999];  % Upper bound on the estimates
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

