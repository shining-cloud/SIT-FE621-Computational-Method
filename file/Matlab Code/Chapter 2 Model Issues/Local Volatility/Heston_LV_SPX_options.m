% Heston Local Volatility on SPX options.
% 1. Estimate one set of parameters for all maturities.
% 2. Build IV and LV for all maturities using this single parameter set.

clc; clear;

% Gauss Laguerre abscissas and weights
[x w] = GenerateGaussLaguerre(32);

% Increments for finite differences for derivatives
dt = 1e-5;
dK = 1e-2;

%% Input the market data
% Spot price, risk free rate, dividend yield
S  = 1164.97;
rf = 0;
q  = 0;

% Maturities and Strikes
T = [0.0192 0.0411 0.1178 0.1945];
K = [1100:5:1230]';
Days = [7 15 43 71];

% Market Implied Volatilities
IV = [...
    0.4232	0.3817	0.3779	0.3673;    0.4195	0.3782	0.3751	0.3644;
    0.4150	0.3751	0.3719	0.3619;    0.4118	0.3715	0.3688	0.3601;
    0.4076	0.3682	0.3655	0.3576;    0.4052	0.3656	0.3630	0.3549;
    0.4007	0.3631	0.3597	0.3524;    0.3967	0.3600	0.3579	0.3499;
    0.3940	0.3557	0.3541	0.3477;    0.3898	0.3526	0.3516	0.3459;
    0.3850	0.3476	0.3470	0.3425;    0.3825	0.3459	0.3455	0.3405;
    0.3777	0.3422	0.3438	0.3384;    0.3750	0.3411	0.3396	0.3354;
    0.3714	0.3382	0.3370	0.3335;    0.3690	0.3349	0.3344	0.3312;
    0.3640	0.3347	0.3315	0.3288;    0.3594	0.3297	0.3292	0.3264;
    0.3570	0.3266	0.3260	0.3246;    0.3552	0.3243	0.3256	0.3221;
    0.3533	0.3216	0.3211	0.3210;    0.3524	0.3195	0.3210	0.3181;
    0.3518	0.3176	0.3188	0.3164;    0.3539	0.3142	0.3167	0.3157;
    0.3556	0.3124	0.3144	0.3139;    0.3615	0.3120	0.3103	0.3121;
    0.3745	0.3120	0.3097	0.3089];

%% Heston parameters on all volatilities combined (one set of parameters)
param = [0.12582  0.001998   0.02272   0.7741  -0.9555];
v0    = param(1); 
theta = param(2); 
kappa = param(3); 
sigma = param(4); 
rho   = param(5);
lambda = 0;
 
% Obtain the local volatilies using a single set of parameters,
% and using finite differences
trap=1;
for t=1:length(T)
	for k=1:length(K)
		LVFD(k,t) = HestonLVFD(S,K(k),T(t),rf,kappa,theta,sigma,lambda,v0,rho,trap,x,w,dt,dK);
		LVAN(k,t) = HestonLVAnalytic(S,K(k),T(t),rf,kappa,theta,sigma,lambda,v0,rho,x,w,trap);
		LVAP(k,t) = HestonLVApprox(S,K(k),T(t),kappa,theta,sigma,v0,rho);
	end
end


%% Plot the Market IV, Heston IV, Heston LV, single parameter set
for t=1:length(T)
	subplot(2,2,t);
	plot(K,LVFD(:,t),'kx',...
		 K,LVAN(:,t),'ko',...
		 K,LVAP(:,t),'k-',...
		 K,IV(:,t)  ,'rx-')
	legend('LV F.D.', 'LV Analytic', 'LV Approx', 'Market IV')
    if t==1
        legend('LV F.D.', 'LV Analytic', 'LV Approx', 'Market IV','Location','SouthWest')
    end
	title(['Maturity ' num2str(Days(t)) ' days']);
end


