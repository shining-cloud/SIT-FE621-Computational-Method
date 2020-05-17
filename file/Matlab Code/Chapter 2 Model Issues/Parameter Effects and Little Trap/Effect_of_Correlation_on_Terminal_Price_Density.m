% Illustration of the log-returns density
% Effect of correlation (rho)

clc; clear;

% Option settings from Table 1 of Heston (1993);
S = 100;
K = 100;            % Strike price
T = 0.5;            % Maturity in years
rf = 0.0;           % Risk free rate
q = 0.0;
x0 = log(S);

% Heston parameters, correlation = -0.5
paramneg.kappa =  2.0;
paramneg.theta =  0.01;
paramneg.sigma =  0.1;
paramneg.v0    =  0.01;
paramneg.rho   = -0.8;
trap = 1;
lambda = 0;

% Heston parameters, correlation = 0
param0.kappa =  2.0;
param0.theta =  0.01;
param0.sigma =  0.1;
param0.v0    =  0.01;
param0.rho   =  0.0;
trap = 1;
lambda = 0;

% Heston parameters, correlation = +0.5
parampos.kappa =  2.0;
parampos.theta =  0.01;
parampos.sigma =  0.1;
parampos.v0    =  0.01;
parampos.rho   =  0.8;
trap = 1;
lambda = 0;


%% Obtain the density of the returns log(S(T)/S(0)) by inverting the characteristic function
dphi = 0.01;
phi = [1e-10:dphi:100];
dx = 0.005;
xT = [4.3:dx:4.9];
for x=1:length(xT)
    intneg = real(exp(-i.*phi.*xT(x)).*HestonCF(phi,paramneg,T,S,rf,q,trap));
    int0   = real(exp(-i.*phi.*xT(x)).*HestonCF(phi,param0  ,T,S,rf,q,trap));
    intpos = real(exp(-i.*phi.*xT(x)).*HestonCF(phi,parampos,T,S,rf,q,trap));
    fneg(x) = trapz(intneg)/pi*dphi;
    f0(x)   = trapz(int0)  /pi*dphi;
    fpos(x) = trapz(intpos)/pi*dphi;
end


%% Plot the densities
plot(xT,fneg,'k--',xT,f0,'r-',xT,fpos,'k+-')
legend('Rho = -0.8', 'Rho = 0.0','Rho = +0.8')
xlabel('Terminal Log Stock Price');
axis([4.3 4.9 0 6])





