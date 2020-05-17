clc; clear;
clear;

% Spot price, risk free rate, dividend yield
S  = 100;
K  = 100;
rf = 0.0;
q  = 0.0;
trap = 1;
kappa = 3;
lambda = 0;
v0 = 0.04;
PutCall = 'C';

% Construct the parameters
for i=1:40;
	T(i)     =  i/4;
  	theta(i) =  0.04 + (i-1)*0.05/100;
  	sigma(i) =  0.30 + (i-1)*0.50/100;
  	rho(i)   = -0.20 + (i-1)*0.35/100;
end;
N = length(T);

% Compute maturity intervals for the Mikhailov Nogel model
tau(1) = T(1);
for t=2:N
	tau(t) = T(t) - T(t-1);
end

% Gauss Laguerre abscissas and weights
[x,w] = GenerateGaussLaguerre(32);


%% Piecewise constant models
% BGM (2010)
param = [kappa v0];
for t=1:40;
    param = [param theta(t) sigma(t) rho(t)];
	ApproxPW(t) = BGMApproxPriceTD(param,T(1:t),S,K,rf,q,PutCall);
end

% True BGM (2010) prices from their paper
TruePW = [3.93 5.53 7.85 11.23 13.92 18.37 22.15 27.17];

% Mikhailov-Nogel (2003)
MNparam0 = [];
tau0 = [];
for t=1:40;
	if t==1
		ClosedPW(t) = HestonPriceGaussLaguerre(PutCall,S,K,T(1),rf,q,kappa,theta(1),sigma(1),lambda,v0,rho(1),trap,x,w);
	else
		MNparam0 = [kappa theta(t-1) sigma(t-1) v0 rho(t-1); MNparam0];
		tau0 = [tau(t-1) tau0];
		MNparam = [kappa theta(t) sigma(t) v0 rho(t)];
		ClosedPW(t) = HestonPriceGLTD(MNparam,MNparam0,tau(t),tau0,K,S,PutCall,rf,q,x,w);
	end
end

% Select the put prices only at the desired maturities
I = ismember(T,[3/12 6/12 1 2 3 5 7 10]);
ApproxPW = ApproxPW(I);
ClosedPW = ClosedPW(I);

%% Heston closed form model with averaged parameters
clear theta sigma rho
T = [3/12 6/12 1 2 3 5 7 10];
N = length(T);

% The averaged parameter values.  Note: v0(5) has been changed
v0     =  [.04 .0397 .0328 .0464 .05624 .2858 .8492 .1454]';
theta  =  [.04 .0404 .0438 .0402 .0404  .0268 .0059 .0457]';
sigma  =  [.30 .3012 .3089 .3112 .3210  .3363 .3541 .3998]';
rho    = -[.20 .1993 .1972 .1895 .1820  .1652 .1480 .1232]';

% The time-static prices using the averaged parameter values
[x,w] =  GenerateGaussLaguerre(32);
for t=1:N;
	ClosedAvg(t) = HestonPriceGaussLaguerre(PutCall,S,K,T(t),rf,q,kappa,theta(t),sigma(t),lambda,v0(t),rho(t),trap,x,w);
end
error = TruePW - ApproxPW;
abserror = sum(abs(TruePW - ApproxPW));

%% Display the results
fprintf('Prices of puts in the Heston model with time-dependent parameters\n')
fprintf('Reproduces column 4 (ATM only) of Table 8 of Benhamou et al. (2010)\n');
fprintf('----------------------------------------------------------------------------------\n');
fprintf('           Heston Closed | Mikhailov-Nogel |  BGM Matlab  | BGM True Price |\n')
fprintf('Maturity  Averaged Param |   Piecewise     |   Piecewise  |   Piecewise    | error\n')
fprintf('----------------------------------------------------------------------------------\n');
for t=1:N
    fprintf('%5.2f %12.2f %18.2f %18.2f %14.2f %10.2f\n',T(t),ClosedAvg(t),ClosedPW(t),ApproxPW(t),TruePW(t),error(t));
end
fprintf('----------------------------------------------------------------------------------\n');
fprintf('Sum absolute error %10.8f \n', abserror)
fprintf('----------------------------------------------------------------------------------\n');


