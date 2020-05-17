% Dougle Heston parameter estimation using loss functions.
% Uses DIA options on May 10, 2012

clc; clear;

% Black Scholes call and put
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));
  
  
%% Load the DIA data
PutIV = [...
    19.62	19.47	20.19	21.15;
    19.10	19.05	19.8	20.82;
    18.60	18.61	19.43	20.57;
    18.10	18.12	19.07	20.21;
    17.61	17.64	18.71	20.00;
    17.18	17.43	18.42	19.74;
    16.71	17.06	18.13	19.50;
    16.44	16.71	17.83	19.27;
    16.61	16.41	17.60	18.99;
    16.61	16.25	17.43	18.84;
    17.01	16.02	17.26	18.71;
    17.55	16.10	17.16	18.46;
    17.86	16.57	17.24	18.42]./100;
K = (124:136);
T = [37 72 135 226]./365;
MktIV = PutIV;
[NK NT] = size(MktIV);
S = 129.14;
rf = 0.0010;
q  = 0.0068;
PutCall = repmat('P',NK,NT);


%% Find the market prices
MktPrice = zeros(NK,NT);
for k=1:NK
    for t=1:NT
        if PutCall(k,t) == 'C';
            MktPrice(k,t) = BSC(S,K(k),rf,q,MktIV(k,t),T(t));
        else 
            MktPrice(k,t) = BSP(S,K(k),rf,q,MktIV(k,t),T(t));
        end
    end
end

%% Settimgs for the Parameter Estimation
% Initial values for kappa,theta,sigma,v0,rho
start = [9 .1 .1 .1 -.9];

% Select the objective function
ObjFun = 1;

% Gauss Laguerre Weights and abscissas
[x w] = GenerateGaussLaguerre(32);
trap = 1;

% Bounds on the parameter estimates
e = 1e-5;
lb = [e   e  e  e -.999];  % Lower bound on the estimates
ub = [20  2  2  3  .999];  % Upper bound on the estimates

%% Single Heston parameter estimates
tic
[param feval] = fmincon(@(p) HestonObjFun(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap),start,[],[],[],[],lb,ub);
t1 = toc;

% Parameter estimates
kappa  = param(1);
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);
lambda = 0;


%% Double Heston parameter estimates
start2 = [start start];
lb2 = [lb lb];      % Lower bound on the estimates
ub2 = [ub ub];      % Upper bound on the estimates

% Ordinary objective function
tic
[Dparam f] = fmincon(@(p) DoubleHestonObjFun(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap),start2,[],[],[],[],lb2,ub2);
t2 = toc;

% Strike Vector Computation objective function
tic
[DparamS fS] = fmincon(@(p) DoubleHestonObjFunSVC(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,ObjFun,x,w,trap),start2,[],[],[],[],lb2,ub2);
t3 = toc;


%% Fit the model implied volatilities from the model prices
% Settings for the bisection method
a = .01;
b = 2;
Tol = 1e-5;
MaxIter = 5000;

for k=1:NK
	for t=1:NT
 		SinglePrice(k,t) = HestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,param,trap,x,w);
		DoublePrice(k,t) = DoubleHestonPriceGaussLaguerre(PutCall(k,t),S,K(k),T(t),rf,q,Dparam,x,w,trap);
		SingleIV(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,SinglePrice(k,t),Tol,MaxIter);
		DoubleIV(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,DoublePrice(k,t),Tol,MaxIter);
		SError(k,t)  = (MktIV(k,t) - SingleIV(k,t))^2;
		DError(k,t)  = (MktIV(k,t) - DoubleIV(k,t))^2;
	end
end

%% Estimation Errors
SingleError = sqrt(sum(sum(SError))) /(NK*NT);
DoubleError = sqrt(sum(sum(DError))) /(NK*NT);

%% Display the results
clc;
fprintf('Method          kappa      theta      sigma        v0        rho   Est Time   IVMSE\n')
fprintf('======================================================================================\n')
fprintf('Univariate  %10.4f %10.4f %10.4f %10.4f %10.4f %6.2f %12.3e\n',param(1), param(2), param(3), param(4),  param(5),t1,SingleError)
fprintf('--------------------------------------------------------------------------------------\n')
fprintf('Double      %10.4f %10.4f %10.4f %10.4f %10.4f %6.2f %12.3e\n',Dparam(1),Dparam(2),Dparam(3),Dparam(4),Dparam(5),t2,DoubleError)
fprintf('Ordinary    %10.4f %10.4f %10.4f %10.4f %10.4f      \n',Dparam(6),Dparam(7),Dparam(8),Dparam(9),Dparam(10))
fprintf('--------------------------------------------------------------------------------------\n')
fprintf('Double      %10.4f %10.4f %10.4f %10.4f %10.4f %6.2f\n',DparamS(1),DparamS(2),DparamS(3),DparamS(4),DparamS(5),t3)
fprintf('SVC ObjFun  %10.4f %10.4f %10.4f %10.4f %10.4f      \n',DparamS(6),DparamS(7),DparamS(8),DparamS(9),DparamS(10))
fprintf('--------------------------------------------------------------------------------------\n')


%% Plot the implied volatilities
SingleIV(find(SingleIV<0)) = NaN;
DoubleIV(find(DoubleIV<0)) = NaN;

for t=1:NT
	subplot(2,2,t)
 	plot(K,MktIV(:,t),'ko-',K,SingleIV(:,t),'rx-',K,DoubleIV(:,t),'bx-')
	title(['Maturity ' num2str(T(t)*365) ' days'])
	legend('Market', 'Single', 'Double')
    if t==NT
        legend('Market', 'Single', 'Double','Location','SouthWest')
    end
    axis([124 136 0.15 0.22])
end

%% plot the implied volatility and local volatility surface
Mat = [T(1):.025:T(end)];
Strike = [K(1):.5:K(end)];
dt = 0.001;
dK = 0.01;
for t=1:length(Mat);
    for k=1:length(Strike);
		Price = DoubleHestonPriceGaussLaguerre('C',S,Strike(k),Mat(t),rf,q,Dparam,x,w,trap);
		IV(k,t) = BisecBSIV('C',S,Strike(k),rf,q,Mat(t),a,b,Price,Tol,MaxIter);
		CT  = DoubleHestonPriceGaussLaguerre('C',S,Strike(k),Mat(t)+dt,rf,q,Dparam,x,w,trap);
		CT_ = DoubleHestonPriceGaussLaguerre('C',S,Strike(k),Mat(t)-dt,rf,q,Dparam,x,w,trap);
        dCdT = (CT - CT_) / (2*dt);
		CK  = DoubleHestonPriceGaussLaguerre('C',S,Strike(k)+dK,Mat(t),rf,q,Dparam,x,w,trap);
		CK0 = DoubleHestonPriceGaussLaguerre('C',S,Strike(k)   ,Mat(t),rf,q,Dparam,x,w,trap);
		CK_ = DoubleHestonPriceGaussLaguerre('C',S,Strike(k)-dK,Mat(t),rf,q,Dparam,x,w,trap);
        dC2dK2 = (CK - 2*CK0 + CK_) / (dK)^2;
        LocalVar = 2*dCdT / (Strike(k)^2*dC2dK2);
        % Local volatility
        LV(k,t) = sqrt(LocalVar);
    end
end

%% Transparent mesh for implied volatility
subplot(1,1,1);
mesh(IV)
xlabel('Maturity')        
ylabel('Strike Price')
set(gca,'XTick',      [0    5    10   15   20   25]);      % Maturity
set(gca,'XTickLabel', [0.1  0.2  0.3  0.4  0.5  0.6])        
set(gca,'YTick',      [0   5   10  15  20  25  30]); % Strike
set(gca,'YTickLabel', [124 126 128 130 132 134 136]);
alpha(0.5);
zlim([0.1 .25]);

% Surface plot for local volatility
hold on
surf(LV)
hold off
