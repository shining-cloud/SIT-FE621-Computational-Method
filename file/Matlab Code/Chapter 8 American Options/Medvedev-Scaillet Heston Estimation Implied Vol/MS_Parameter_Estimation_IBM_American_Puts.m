% Heston parameter estimation using implied vols of American put prices on IBM
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

% Option features
T = days./365;
S = 122.10;
rf = 0.015;
q  = 0.01;


%% Select a single maturity only
t = 1;
MktPrice = MktPrice(:,t);
T = T(t);
days = days(t);

% Dimensions of the prices
[NK NT] = size(MktPrice);
PutCall = 'P';


%% Implied vol by MS Black Scholes put and bisection algorithm
a = 0.01;
b = 2;
Tol = 1e-5;
MaxIter = 1000;

% American implied volatilities
for t=1:NT
    for k=1:NK
        MktIV(k,t) = BisecMSIV(S,K(k),rf,q,T(t),a,b,MktPrice(k,t),Tol,MaxIter);
    end
end


%% Settings for the M-S put
yinf = 1e4;
NumTerms = 4;

% Settings for the Newton-Cotes method in the M-S put price
method = 3;
A = 1e-20;
B = 100;
N = 1000;
trap = 1;


%% Parameter estimation
% kappa theta sigma v0 rho
start = [20.0000, 0.0357, 3.9103, 0.1586, -0.4611];

% Estimation bounds
e = 1e-3;
lb = [e   e  e  e -.99];
ub = [20  2  5  2  .99];

% Bisection algorithm settings
a = 0.01;
b = 2.0;
Tol = 1e-8;
MaxIter = 2000;

% Run the Nelder-Mead algorithm with parameter bounds
fprintf('------------------------------------------------------------------\n')
fprintf('LossFn          kappa      theta       sigma       v0        rho  \n')
fprintf('------------------------------------------------------------------\n')
[param feval] = fmincon(@(p) HestonObjFunMSIV(p,S,rf,q,K,T,trap,MktIV,method,A,B,N,yinf,NumTerms,a,b,Tol,MaxIter),start,[],[],[],[],lb,ub);
fprintf('------------------------------------------------------------------\n')


%% Find MS American implied volatilities values with the estimated parameters
param = [18.5004, 0.03534, 3.93462, 0.16027, -0.50169];
for k=1:NK
	for t=1:NT
         ModelPrice(k,t) = MSPrice(S,K(k),T(t),rf,q,param,trap,method,A,B,N,NumTerms,yinf);
         ModelIV(k,t) = BisecMSIV(S,K(k),rf,q,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
         error(k,t) = (ModelIV(k,t) - MktIV(k,t))^2 ;
	end
end
MSE = sum(sum(error)) / (NT*NK);

%% Display the results
fprintf('  \n')
fprintf([num2str(NumTerms) '-term Heston-Medvedev-Scaillet parameter estimation \n'])
fprintf('--------------------------------------------------------\n')
fprintf('  kappa    theta    sigma      v0       rho     MSE \n')
fprintf('--------------------------------------------------------\n')
fprintf('%8.4f %8.4f %8.4f %8.4f %8.4f %12.4e \n',param(1),param(2),param(3),param(4),param(5),MSE) 
fprintf('--------------------------------------------------------\n')

%% Plot the results
plot(K,MktIV,'ko-',K,ModelIV,'ro-')
xlabel('Strike Price')
ylabel('Implied Volatility')
legend('Market IV','Model IV')

