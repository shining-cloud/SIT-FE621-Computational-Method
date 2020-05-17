% Illustration of the Elices Time dependent Heston model

clc; clear;

% Black Scholes call and put
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));
  
%% Quoted prices and implied vols
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
PutCall = 'P';

% Maturity increments
tau(1) = T(1);
for t=2:NT;
    tau(t) = T(t) - T(t-1);
end

% Weights and abscissas
[x w] = GenerateGaussLaguerre(32);
trap = 1;

%% Find the market prices
for k=1:NK
    for t=1:NT
        if PutCall == 'C';
            MktPrice(k,t) = BSC(S,K(k),rf,q,MktIV(k,t),T(t));
        else 
            MktPrice(k,t) = BSP(S,K(k),rf,q,MktIV(k,t),T(t));
        end
    end
end
a = 0.0001;
b = 3;
Tol = 1e-10;
MaxIter = 10000;

%% Static parameter estimates, prices, and implied vol
% kappa theta sigma v0 rho
start = [2 0.1 1.2 0.05 -.5];
e = 1e-3;
lb = [e   e   e   e  -.999];  % Lower bound on the estimates
ub = [20  10  10  10  .999];  % Upper bound on the estimates

% Select the loss function
ObjFun = 1;

% Parameter estimation
[ParamTI fe] = fmincon(@(b) HestonObjFun(b,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap), start,[],[],[],[],lb,ub);

% Prices and implied vol
for t=1:NT
    for k=1:NK
        PriceTI(k,t) = HestonPriceGaussLaguerre(PutCall,S,K(k),T(t),rf,q,ParamTI,trap,x,w);
        IVTI(k,t) = BisecBSIV(PutCall,S,K(k),rf,q,T(t),a,b,PriceTI(k,t),Tol,MaxIter);
    end
    ErrorTI(t) = sum((IVTI(:,t) - MktIV(:,t)).^2) / NK;
end

% Volatility estimated from first maturity only
[ParamTemp fe] = fmincon(@(b) HestonObjFun(b,S,rf,q,MktPrice(:,1),K,T(1),PutCall,MktIV(:,1),ObjFun,x,w,trap), start,[],[],[],[],lb,ub);
v0 = ParamTemp(4);


%% Time dependent parameter estimates, prices, and implied vol
% Remove v0 from subsequent parameter estimation,
% update lower and upper parameter bounds
start(4) = [];
lb(4)    = [];
ub(4)    = [];

% Elices (2009) model for remaining parameters
[ParamTD(1,:) fe] = fmincon(@(b) ElicesObjFun(b,[]          ,  v0,tau(1)  ,S,K,rf,q,MktPrice(:,1),MktIV(:,1),PutCall,ObjFun,x,w,trap), start,[],[],[],[],lb,ub);
start = ParamTD(1,:);
[ParamTD(2,:) fe] = fmincon(@(b) ElicesObjFun(b,ParamTD(1,:),  v0,tau(1:2),S,K,rf,q,MktPrice(:,2),MktIV(:,2),PutCall,ObjFun,x,w,trap), start,[],[],[],[],lb,ub);
start = ParamTD(2,:);
[ParamTD(3,:) fe] = fmincon(@(b) ElicesObjFun(b,ParamTD(1:2,:),v0,tau(1:3),S,K,rf,q,MktPrice(:,3),MktIV(:,3),PutCall,ObjFun,x,w,trap), start,[],[],[],[],lb,ub);
start = ParamTD(3,:);
[ParamTD(4,:) fe] = fmincon(@(b) ElicesObjFun(b,ParamTD(1:3,:),v0,tau(1:4),S,K,rf,q,MktPrice(:,4),MktIV(:,4),PutCall,ObjFun,x,w,trap), start,[],[],[],[],lb,ub);


%% Prices and implied vol
clc;
 for k=1:NK
     PriceTD(k,1) = ElicesPrice(PutCall,S,K(k),tau(1)  ,rf,q,ParamTD(1,:),[]            ,v0,trap,x,w);
     PriceTD(k,2) = ElicesPrice(PutCall,S,K(k),tau(1:2),rf,q,ParamTD(2,:),ParamTD(1,:)  ,v0,trap,x,w);
     PriceTD(k,3) = ElicesPrice(PutCall,S,K(k),tau(1:3),rf,q,ParamTD(3,:),ParamTD(1:2,:),v0,trap,x,w);
     PriceTD(k,4) = ElicesPrice(PutCall,S,K(k),tau(1:4),rf,q,ParamTD(4,:),ParamTD(1:3,:),v0,trap,x,w);
end
for t=1:NT
    for k=1:NK
        IVTD(k,t) = BisecBSIV(PutCall,S,K(k),rf,q,T(t),a,b,PriceTD(k,t),Tol,MaxIter);
    end
    ErrorTD(t) = sum((IVTD(:,t) - MktIV(:,t)).^2) / NK;
end

% IVRMSE errors
SumErrorTI = sum(ErrorTI)/NT;
SumErrorTD = sum(ErrorTD)/NT;


%% Display the results
fprintf('Elices Time Dependent Parameter Estimates with loss function %1.0f\n',ObjFun)
fprintf('Maturity   kappa      theta      sigma       rho        v0      IVMSE\n')
fprintf('-------------------------------------------------------------------------\n');
for t=1:NT
    fprintf('%6.0f %10.4f %10.4f %10.4f %10.4f %10.4f %10.2d\n',T(t)*365,ParamTD(t,:),v0,ErrorTD(t));
end
fprintf('-------------------------------------------------------------------------\n');
fprintf('Static Parameter Estimates\n')
fprintf('Maturity   kappa      theta      sigma       v0        rho\n')
fprintf('-------------------------------------------------------------\n');
fprintf('  All  %10.4f %10.4f %10.4f %10.4f %10.4f \n',ParamTI);
fprintf('-------------------------------------------------------------\n');
fprintf('IVMSE from time dependent model  %10.2d \n',SumErrorTD);
fprintf('IVMSE from time static model     %10.2d \n',SumErrorTI);
fprintf('-------------------------------------------------------------\n');


%% Plot the implied volatilities
for t=1:NT
	subplot(2,2,t)
 	plot(K,MktIV(:,t),'ko-',K,IVTD(:,t),'rx-',K,IVTI(:,t),'bx-')
	legend('Market', 'Elices','Static')
	title(['Maturity ' num2str(T(t)*365) ' days']); 
    axis([124 136 0.15 0.22])
end
