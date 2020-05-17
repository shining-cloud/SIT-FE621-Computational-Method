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
T = [37 72 135 226 278]./365;
S = 129.14;
rf = 0.0010;
q  = 0.0068;

MktIV = PutIV;
[NK NT] = size(MktIV);
PutCall = repmat('P',NK,NT);

%% Find the market prices
for k=1:NK
    for t=1:NT
        if PutCall(k,t) == 'C';
            MktPrice(k,t) = BSC(S,K(k),rf,q,MktIV(k,t),T(t));
        else 
            MktPrice(k,t) = BSP(S,K(k),rf,q,MktIV(k,t),T(t));
        end
    end
end

% Weights and abscissas
[x w] = GenerateGaussLaguerre(32);

%% Create the maturity increments
tau(1) = T(1);
for t=2:NT
	tau(t) = T(t) - T(t-1);
end

%% Combined estimates for Time Independent parameters -- starting values
% parameter order [kappa theta sigma v0 rho]
start = [2 .1 1.2 .05 -.5];
e = 1e-3;
lb = [e   e   e   e  -.999];  % Lower bound on the estimates
ub = [20  10  10  10  .999];  % Upper bound on the estimates

trap = 1;
ObjFun = 3;
[ParamTI fe] = fmincon(@(b) HestonObjFun(b,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap),start,[],[],[],[],lb,ub)


%% Time Dependent parameters
% Initialize the matrix of time dependent parameters
ParamTD = zeros(NT,5);

% First set of time-dependent parameters (t=1)
[ParamTD f] = fmincon(@(param) MNObjFun(param,[],tau(1),[],S,K,rf,q,MktPrice(:,1),PutCall(:,1),MktIV(:,1),ObjFun,x,w), start,[],[],[],[],lb,ub)
oldparam = [];
oldtau   = [];

%% Remaining sets of time-dependent parameters (t=2 to t=T)
% parameter order [kappa(1) theta(1) sigma(1) v0(1) rho(1); 
%                  kappa(2) theta(2) sigma(2) v0(2) rho(2); etc.]
for t=2:NT
	oldparam = [oldparam; ParamTD(t-1,:)];
	oldtau   = [oldtau, tau(t-1)];
	start = ParamTD(t-1,:);
	[ParamTD(t,:) f] = fmincon(@(param) MNObjFun(param,oldparam,tau(t),oldtau,S,K,rf,q,MktPrice(:,t),PutCall(:,t),MktIV(:,t),ObjFun,x,w), start,[],[],[],[],lb,ub)
end

%%  Fit Heston prices to both sets of parameters
% Settings for the bisection method
Tol = 1e-7;
MaxIter = 1000;
a = 0.001;
b = 5;

for k=1:NK
	% Initialize past parameter estimates and past maturities
	oldparam = [];
	oldtau   = [];
	for t=1:NT
		% Prices using time-independent parameters floored at $0.01
		PriceTI(k,t) = MNPriceGaussLaguerre(ParamTI,[],T(t),[],K(k),S,PutCall(k,t),rf,q,x,w);
		PriceTI(k,t) = max(0.01,PriceTI(k,t));
		% Prices using time-dependent parameters floored at $0.01
		if t==1
			PriceTD(k,t) = MNPriceGaussLaguerre(ParamTD(t,:),[],tau(t),[],K(k),S,PutCall(k,t),rf,q,x,w);
		else
			oldparam = [oldparam; ParamTD(t-1,:)];
			oldtau   = [oldtau, tau(t-1)];
			PriceTD(k,t) = MNPriceGaussLaguerre(ParamTD(t,:),oldparam,tau(t),oldtau,K(k),S,PutCall(k,t),rf,q,x,w);
		end
		PriceTD(k,t) = max(0.01,PriceTD(k,t));
		% Implied volatilities
		IVTI(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,PriceTI(k,t),Tol,MaxIter);
		IVTD(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,PriceTD(k,t),Tol,MaxIter);
	end
end

ErrorTD = sum(sum((IVTD - MktIV).^2)) / (NT*NK);
ErrorTI = sum(sum((IVTI - MktIV).^2)) / (NT*NK);


%% Display the results
clc;
fprintf('Mikhail Nogel Time Dependent Parameter Estimates\n')
fprintf('Maturity   kappa      theta      sigma      v0        rho\n')
fprintf('-------------------------------------------------------------\n');
for t=1:NT
    fprintf('%6.0f %10.4f %10.4f %10.4f %10.4f %10.4f \n',T(t)*365,ParamTD(t,:));
end
fprintf('-------------------------------------------------------------\n');
fprintf('Constant Parameter Estimates\n')
fprintf('Maturity   kappa      theta      sigma      v0        rho\n')
fprintf('-------------------------------------------------------------\n');
fprintf('  All %10.4f %10.4f %10.4f %10.4f %10.4f \n',ParamTI);
fprintf('-------------------------------------------------------------\n');
fprintf('IVMSE from time dependent model  %10.2d \n',ErrorTD);
fprintf('IVMSE from static model          %10.2d \n',ErrorTI);


%% Plots of implied volatility
for t=1:NT
	subplot(2,2,t)
	plot(K,IVTI(:,t),'bx-',K,IVTD(:,t),'rx-',K,MktIV(:,t),'ko-')
    axis([124 136 0.15 0.22])
	legend('Static', 'M-N', 'Market')
	title(['Maturity ' num2str(T(t)*365) ' days']); 
end

%% Total volatility -- static model'
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

   