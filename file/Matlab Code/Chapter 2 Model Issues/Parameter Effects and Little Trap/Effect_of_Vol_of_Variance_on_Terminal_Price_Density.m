% Illustration of the log-returns density
% Effect of vol of variance (sigma)

clc; clear;

% Option settings from Table 1 of Heston (1993);
S = 100;
K = 100;            % Strike price
T = 0.5;            % Maturity in years
rf = 0.0;           % Risk free rate
q = 0.0;

% Heston parameters, sigma = 0
paramlo.kappa =  2.0;
paramlo.theta =  0.01;
paramlo.sigma =  1e-5;
paramlo.v0    =  0.01;
paramlo.rho   =  0.0;
trap = 1;
lambda = 0;

% Heston parameters, sigma = 0.2
parammed.kappa =  2.0;
parammed.theta =  0.01;
parammed.sigma =  0.20;
parammed.v0    =  0.01;
parammed.rho   =  0.0;
trap = 1;
lambda = 0;

% Heston parameters, sigma = 0.4
paramhi.kappa =  2.0;
paramhi.theta =  0.01;
paramhi.sigma =  0.40;
paramhi.v0    =  0.01;
paramhi.rho   =  0.0;
trap = 1;
lambda = 0;

%% Obtain the density of log(S(T)) by inverting the characteristic function
dphi = 0.01;
phi = [1e-10:dphi:100];
dx = 0.001;
xT = [4.3:dx:4.9];
for x=1:length(xT)
    intlo  = real(exp(-i.*phi.*xT(x)).*HestonCF(phi,paramlo ,T,S,rf,q,trap));
    intmed = real(exp(-i.*phi.*xT(x)).*HestonCF(phi,parammed,T,S,rf,q,trap));
    inthi  = real(exp(-i.*phi.*xT(x)).*HestonCF(phi,paramhi ,T,S,rf,q,trap));
    flo(x)  = trapz(intlo) /pi*dphi;
    fmed(x) = trapz(intmed)/pi*dphi;
    fhi(x)  = trapz(inthi) /pi*dphi;
end

%% Plot the densities
plot(xT,flo,'k--',xT,fmed,'r-',xT,fhi,'k-')
legend('Sigma = 0.0', 'Sigma = 0.2','Sigma = 0.4')
xlabel('Terminal Log Stock Price');
axis([4.3 4.9 0 9])


