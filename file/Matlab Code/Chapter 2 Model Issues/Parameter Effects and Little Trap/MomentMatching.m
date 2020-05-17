% Matched moment for the Heston model to retrieve BS volatility

clc; clear;

% Parameter values from Table 1 of Heston (1993)
param.kappa = 2;
param.theta = 0.01;
param.sigma = 0.1;
param.v0 = 0.01;
lambda = 0;

% Maturity, strike, spot, etc.
tau = 0.5;
K = 100;
S = 100;
r = 0;
q = 0;
PutCall = 'C';
trap = 1;

% Different values of rho
rho = [-0.5 0.0 0.5];

% Find the moments of the Heston c.f. by finite differences
dphi = 1e-4;
fprintf('      Rho        BSvol     Matched BSvol\n')
fprintf('----------------------------------------\n')
for k=1:3
    param.rho = rho(k);
    dfp = HestonCF(+dphi,param,tau,S,r,q,trap);
    dfm = HestonCF(-dphi,param,tau,S,r,q,trap);
    df = (dfp - dfm)/2/dphi;
    EX = df/i;
    dfpp = HestonCF(+2*dphi,param,tau,S,r,q,trap);
    dfmm = HestonCF(-2*dphi,param,tau,S,r,q,trap);
    dff1 = (dfpp - dfp)/dphi;
    dff2 = (dfm - dfmm)/dphi;
    ddf  = (dff1 - dff2)/3/dphi;
    EX2 = ddf/i^2;
    var = EX2 - EX^2;
    BSvol(k) = sqrt(var);
    fprintf('%10.2f   %10.4f   %10.4f \n',rho(k),BSvol(k),BSvol(k).*sqrt(2))
end;


