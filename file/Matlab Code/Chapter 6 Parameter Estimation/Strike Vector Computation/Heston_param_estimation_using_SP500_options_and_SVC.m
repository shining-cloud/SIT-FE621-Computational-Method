% Parameter estimation using loss functions and Strike Vector Computation
% Uses SP500 options

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

%% Parameter Estimation
% Initial values for kappa,theta,sigma,v0,rho
start = [9.0 0.05 0.3 0.05 -0.8];

% Specify the objective function and method to obtain the option price
ObjFun = 1;
Method = 1;

% Gauss Laguerre Weights and abscissas
[x w] = GenerateGaussLaguerre(32);

% Settings for the parameter estimates
a = .01;
b = 2;
Tol = 1e-5;
MaxIter = 1000;
trap = 1;

e = 1e-5;
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [20  2  2  3  .999];  % Upper bound on the estimates


%% Parameter estimates using ordinary objective function
tic
paramORD = fmincon(@(p) HestonObjFun(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,Method,a,b,Tol,MaxIter),start,[],[],[],[],lb,ub);
t1 = toc;

%% Parameter estimates using SVC objective function and Heston c.f.
tic
type = 'Heston';
paramSVC1 = fmincon(@(p) HestonObjFunSVC(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,type),start,[],[],[],[],lb,ub);
t2 = toc;

%% Parameter estimates using SVC objective function and Attari c.f.
tic
type = 'Attari';
paramSVC2 = fmincon(@(p) HestonObjFunSVC(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,trap,ObjFun,a,b,Tol,MaxIter,type),start,[],[],[],[],lb,ub);
t3 = toc;

%% Display the results
clc;
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('Method                     kappa     theta        sigma       v0         rho     Est time\n')
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('Ordinary obj fun       %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',paramORD(1),paramORD(2),paramORD(3),paramORD(4),paramORD(5),t1);
fprintf('SVC obj fun, Heston CF %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',paramSVC1(1),paramSVC1(2),paramSVC1(3),paramSVC1(4),paramSVC1(5),t2);
fprintf('SVC obj fun, Attari CF %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',paramSVC2(1),paramSVC2(2),paramSVC2(3),paramSVC2(4),paramSVC2(5),t3);
fprintf('--------------------------------------------------------------------------------------------\n');

