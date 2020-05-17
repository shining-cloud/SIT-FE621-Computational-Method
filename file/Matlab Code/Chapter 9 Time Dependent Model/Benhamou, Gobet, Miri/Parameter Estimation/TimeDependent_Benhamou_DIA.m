% Benhamou, Gobet, Miri Time-Dependent Heston model
% Parameter estimation using DJIA data
clc; clear;

% Black Scholes call and put
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));
  
%% Load the DIA data
PutIV = [...
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
Days = [37 72 135 226];
T = Days./365;
MktIV = PutIV;
[NK NT] = size(MktIV);
S = 129.14;
rf = 0.0010;
q  = 0.0068;
PutCall = repmat('P',NK,NT);
  

%% Market Prices
for k=1:NK
	for t=1:NT
		if strcmp(PutCall(k,t),'C')
			MktPrice(k,t) = BSC(S,K(k),rf,q,MktIV(k,t),T(t));
		elseif strcmp(PutCall(k,t),'P')
			MktPrice(k,t) = BSP(S,K(k),rf,q,MktIV(k,t),T(t));
		end
	end
end

%% Settings for the objective function
[x w] = GenerateGaussLaguerre(32);
ObjFun = 3;
trap = 1;

%% Combined estimates for Time Independent parameters -- starting values
% Parameter order:[kappa v0 theta sigma  rho]
% Starting values
kappa0 =  2.0;
theta0 =  0.1;
sigma0 =  1.2;
v00    =  0.05;
rho0   = -0.5;

start = [kappa0 theta0 sigma0 v00 rho0];
e = 1e-3;
lb = [e   e  e  e  -.99];  % Lower bound on the estimates
ub = [20  3  3  3   .99];  % Upper bound on the estimates

[ParamTI f] = fmincon(@(p) HestonObjFun(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap), start,[],[],[],[],lb,ub);

%% Time Dependent parameters
% Parameter order:
% [kappa v0 | theta(1) sigma(1) rho(1) | theta(2) sigma(2) rho(2) | etc...]
kappa0 = ParamTI(1);
theta0 = ParamTI(2);
sigma0 = ParamTI(3);
v00    = ParamTI(4);
rho0   = ParamTI(5);

% Starting values
start = [kappa0 v00 repmat([theta0 sigma0 rho0],1,length(T))];

% Parameter bounds
lb = [e   e repmat([e  e  -.999],1,length(T))];
ub = [100 2 repmat([10 10  .999],1,length(T))];

% Find the time dependent parameter estimates
options = optimset('MaxIter',1e5,'MaxFunEvals',1e7);
[paramTD f] = fmincon(@(p) HestonBGMObjFun(p,T,S,K,rf,q,MktPrice,MktIV,ObjFun,PutCall), start,[],[],[],[],lb,ub);

% Number of time dependent parameters
Nparam = 2 + 3*length(T);

kappa = paramTD(1);
v0    = paramTD(2);
for i=3:3:Nparam
	j = floor(i/3);
	theta(j) = paramTD(i);
	sigma(j) = paramTD(i+1);
	rho(j)   = paramTD(i+2);
end

% Stack the estimated parameters in the matrix ParamTD
% Parameters for the shortest maturity are on top, those for the longest
% maturity on the bottom
for t=1:NT
	ParamTD(t,:) = [kappa v0 theta(t) sigma(t) rho(t)];
end

 
%%  Obtain prices for both sets of parameters
% Settings for the bisection method
Tol = 1e-7;
MaxIter = 1000;
a = 0.00001;
b = 5;
trap = 1;

for k=1:NK
	for t=1:NT
		% Prices using time-independent parameters
        PriceTI(k,t) = HestonPriceGaussLaguerre(PutCall,S,K(k),T(t),rf,q,ParamTI,trap,x,w);
		% Prices using time-dependent parameters
		PriceTD(k,t) = BGMApproxPriceTD(paramTD,T(1:t),S,K(k),rf,q,PutCall);
		% Implied volatilities
  		IVTI(k,t) = BisecBSIV(PutCall,S,K(k),rf,q,T(t),a,b,PriceTI(k,t),Tol,MaxIter);
   		IVTD(k,t) = BisecBSIV(PutCall,S,K(k),rf,q,T(t),a,b,PriceTD(k,t),Tol,MaxIter);
	end
end

% IVRMSE
ErrorTI = sum(sum((IVTI - MktIV).^2)) / (NT*NK);
ErrorTD = sum(sum((IVTD - MktIV).^2)) / (NT*NK);


%% Display the results
clc;
fprintf('Estimation on DIA \n')
fprintf('Constant Parameter Estimates\n')
fprintf('Maturity   kappa      theta      sigma        v0       rho\n')
fprintf('------------------------------------------------------------\n');
fprintf('All    %10.4f %10.4f %10.4f %10.4f %10.4f \n',ParamTI(1),ParamTI(2),ParamTI(3),ParamTI(4),ParamTI(5));  
fprintf('------------------------------------------------------------\n');
fprintf(' \n')
fprintf('Benhamou, Gobet, Miri (2010) Time Dependent Parameter Estimates\n')
fprintf('Maturity   kappa        v0       theta      sigma      rho\n')
fprintf('------------------------------------------------------------\n');
for t=1:NT
    fprintf('%6.0f %10.4f %10.4f %10.4f %10.4f %10.4f \n',Days(t),ParamTD(t,1),ParamTD(t,2),ParamTD(t,3),ParamTD(t,4),ParamTD(t,5))
end
fprintf('------------------------------------------------------------\n');
fprintf('IVRMSE from time dependent model  %10.2d \n',ErrorTD);
fprintf('IVRMSE from static model          %10.2d \n',ErrorTI);
fprintf('------------------------------------------------------------\n');


%% Plot the implied volatilities
Dates = {'3 months'; '1 Year'; '18 Months'; '2 Years'};
for t=1:NT
	subplot(2,2,t)
	plot(K,IVTI(:,t),'bx-',K,IVTD(:,t),'rx-',K,MktIV(:,t),'ko-')
	legend('Static', 'BGM', 'Market')
	title(['Maturity ' num2str(T(t)*365) ' days']); 
    axis([124 136 0.15 0.22])
end


%% Total volatility -- static model
% subplot(1,1,1);
% plot(K,IVTI(:,1)*sqrt(T(1)),...
%      K,IVTI(:,2)*sqrt(T(2)),...
%      K,IVTI(:,3)*sqrt(T(3)),...
%      K,IVTI(:,4)*sqrt(T(4)));
% legend(['Maturity ' num2str(T(1)*365) ' days'], ...
%        ['Maturity ' num2str(T(2)*365) ' days'], ...
%        ['Maturity ' num2str(T(3)*365) ' days'], ...
%        ['Maturity ' num2str(T(4)*365) ' days']);

%% Total volatility -- time dependent model
% plot(K,IVTD(:,1)*sqrt(T(1)),...
%      K,IVTD(:,2)*sqrt(T(2)),...
%      K,IVTD(:,3)*sqrt(T(3)),...
%      K,IVTD(:,4)*sqrt(T(4)));
% legend(['Maturity ' num2str(T(1)*365) ' days'], ...
%        ['Maturity ' num2str(T(2)*365) ' days'], ...
%        ['Maturity ' num2str(T(3)*365) ' days'], ...
%        ['Maturity ' num2str(T(4)*365) ' days']);

