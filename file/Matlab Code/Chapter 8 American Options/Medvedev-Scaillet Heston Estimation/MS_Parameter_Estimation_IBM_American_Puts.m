% Heston American put prices on IBM
% Uses Medvedev-Scaillet expansion

clc; clear;

% IBM Put prices for May 7, 2010 closing price;
MktPrice = [...
     0.130	 0.620	 1.275	 2.950	 4.750;
     0.260	 0.955	 1.830	 3.925	 6.075;
     0.485	 1.500	 2.610	 5.200	 7.625;
     0.995	 2.445	 3.775	 6.800	 9.375;
     2.155	 3.900	 5.475	 8.800	11.550;
     4.525	 6.225	 7.775	11.225	13.975;
     8.375	 9.525	10.850	14.125	16.775;
    13.075	13.600	14.575	17.425	19.900;
    17.975	18.175	18.825	21.150	23.325];

% Strikes, maturities, spot price, risk free rate, dividend yield
K = 100:5:140;
days = [14 42 70 161 259];
T = days./365;
S = 122.10;
rf = 0.015;
q  = 0.01;

% Dimensions of the prices
[NK NT] = size(MktPrice);
PutCall = 'P';

% Settings for the M-S put
NumTerms = 3;
yinf = 1e4;

% Settings for the Newton-Cotes method in the M-S put price
method = 3;
A = 1e-20;
B = 100;
N = 1000;
trap = 1;


%% Parameter estimation
start = [0.0146, 0.0284, 0.4010, 0.1066, -0.5072];

% Estimation bounds
e = 1e-3;
lb = [e   e  e  e -.99];
ub = [20  2  5  2  .99];

% optimization settings
options = optimset('MaxIter',25);

% Run the Nelder-Mead algorithm with parameter bounds
fprintf('------------------------------------------------------------------\n')
fprintf('LossFn          kappa      theta       sigma       v0        rho  \n')
fprintf('------------------------------------------------------------------\n')
[param feval] = fminsearch(@(p) HestonObjFunMS(p,S,rf,q,K,T,trap,MktPrice,method,A,B,N,lb,ub,yinf,NumTerms),start,options);
fprintf('------------------------------------------------------------------\n')


%% Find MS American put values with the estimated parameters
for t=1:NT
    for k=1:NK
        AmerPut(k,t) = MSPrice(S,K(k),T(t),rf,q,param,trap,method,A,B,N,NumTerms,yinf);
    end
end

% Mean Square Error
MSE = sum(sum(abs(AmerPut - MktPrice).^2))/(NT*NK);

%% Display the results
fprintf('  \n')
fprintf([num2str(NumTerms) '-term Heston-Medvedev-Scaillet parameter estimation \n'])
fprintf('--------------------------------------------------------\n')
fprintf('  kappa    theta    sigma      v0       rho     MSE \n')
fprintf('--------------------------------------------------------\n')
fprintf('%8.4f %8.4f %8.4f %8.4f %8.4f %12.4e \n',param(1),param(2),param(3),param(4),param(5),MSE) 
fprintf('--------------------------------------------------------\n')

%% Plot the market and model prices
plot(K,MktPrice(:,1),'c-',K,MktPrice(:,2),'g-',K,MktPrice(:,3),'r-',...
     K,MktPrice(:,4),'b-',K,MktPrice(:,5),'k-',...
     K,AmerPut(:,1),'co--',K,AmerPut(:,2),'go--',K,AmerPut(:,3),'ro--',...
     K,AmerPut(:,4),'bo--',K,AmerPut(:,5),'ko--')
xlabel('Strike Price')
ylabel('American Put Price')
legend([' ' num2str(days(1)) '-day maturity'],...
       [' ' num2str(days(2)) '-day maturity'],...
       [' ' num2str(days(3)) '-day maturity'],...
       [    num2str(days(4)) '-day maturity'],...
       [    num2str(days(5)) '-day maturity'],...
       'Location','NorthWest')
pause

   
   