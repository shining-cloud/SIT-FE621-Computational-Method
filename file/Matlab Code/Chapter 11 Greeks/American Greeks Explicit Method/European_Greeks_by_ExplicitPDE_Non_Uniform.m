% Illustration of Greeks by PDE

clc; clear;

% Strike price, risk free rate, dividend yield, and maturity
K = 100;
r = 0.02;
q = 0;
Mat = 0.15;
PutCall = 'C';
EuroAmer = 'E';

% Heston parameters
kappa =  1.5;
theta =  0.04;
sigma =  0.3;
rho   = -0.9;
v0    =  0.05;
lambda = 0;
params = [kappa theta sigma v0 rho lambda];

% Obtain the price by 2-D interpolation
S0 = 101.52;
V0 = 0.05412;

% Exact Greeks for these parameter settings
DeltaClosed =  0.63777;
GammaClosed =  0.03995;
Vega1Closed = 13.55363;
VannaClosed = -0.37411;
VolgaClosed =  8.25269;
ThetaClosed = -12.3258;

% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0;  Smax = 2*K;
Vmin = 0;  Vmax = 0.5;
Tmin = 0;  Tmax = Mat;

% Number of grid points for the stock, volatility, and maturity
nS = 79;        % Stock price
nV = 39;        % Volatility
nT = 3000;      % Maturity

% The maturity time increment and grid
dt = (Tmax-Tmin)/nT;
T = [0:nT].*dt;


%% Pricing Using a Non-Uniform Grid
% The stock price grid
c = K/5;
dz = 1/nS*(asinh((Smax-K)/c) - asinh(-K/c));
for i=1:nS+1;
	z(i) = asinh(-K/c) + (i-1)*dz;
	S(i) = K + c*sinh(z(i));
end

% The volatility grid
d = Vmax/500;
dn = asinh(Vmax/d)/nV;
for j=1:nV+1
	n(j) = (j-1)*dn;
	V(j) = d*sinh(n(j));
end

% Solve the PDE
[U u] = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer);
clear dv ds

%% PDE Price and Greeks by 2-D interpolation
PricePDE = interp2(V,S,U,V0,S0);

% Delta
dS = 1;
D1 = interp2(V,S,U,V0,S0+dS);
D2 = interp2(V,S,U,V0,S0-dS);
DeltaPDE = (D1-D2)/2/dS;

% Gamma
GammaPDE = (D1 - 2*PricePDE + D2)/dS^2;

% Vega #1
dV = 1e-2;
V1 = interp2(V,S,U,V0+dV,S0);
V2 = interp2(V,S,U,V0-dV,S0);
Vega1PDE = (V1-V2)/2/dV*2*sqrt(V0);

% Vanna
C1 = interp2(V,S,U,V0+dV,S0+dS);
C2 = interp2(V,S,U,V0-dV,S0+dS);
C3 = interp2(V,S,U,V0+dV,S0-dS);
C4 = interp2(V,S,U,V0-dV,S0-dS);
VannaPDE = (C1 - C2 - C3 + C4)/4/dV/dS*2*sqrt(V0);
	
% Volga
dC2 = (V1 - 2*PricePDE + V2)/dV^2;
VolgaPDE = 4*sqrt(V0)*(dC2*sqrt(V0) + Vega1PDE/4/V0);

% Theta
T1 = interp2(V,S,U,V0,S0);  % Obtain U(S,v,T)
T2 = interp2(V,S,u,V0,S0);  % Obtain U(S,v,T-dt)
ThetaPDE = -(T1 - T2)/dt;


%% Closed form Price and Greeks by Finite Differences
trap = 1;
[x w] = GenerateGaussLaguerre(32);
PriceClosed = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);

% Delta
D1 = HestonPriceGaussLaguerre(PutCall,S0+dS,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);
D2 = HestonPriceGaussLaguerre(PutCall,S0-dS,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);
DeltaFD = (D1-D2)/2/dS;

% Gamma
GammaFD = (D1 - 2*PriceClosed + D2)/dS^2;

% Vega #1  
V1 = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0+dV,rho,trap,x,w);
V2 = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0-dV,rho,trap,x,w);
Vega1FD = (V1-V2)/2/dV*2*sqrt(V0);

% Vanna
C1 = HestonPriceGaussLaguerre(PutCall,S0+dS,K,Mat,r,q,kappa,theta,sigma,lambda,V0+dV,rho,trap,x,w);
C2 = HestonPriceGaussLaguerre(PutCall,S0+dS,K,Mat,r,q,kappa,theta,sigma,lambda,V0-dV,rho,trap,x,w);
C3 = HestonPriceGaussLaguerre(PutCall,S0-dS,K,Mat,r,q,kappa,theta,sigma,lambda,V0+dV,rho,trap,x,w);
C4 = HestonPriceGaussLaguerre(PutCall,S0-dS,K,Mat,r,q,kappa,theta,sigma,lambda,V0-dV,rho,trap,x,w);
VannaFD = (C1 - C2 - C3 + C4)/4/dV/dS*2*sqrt(V0);
	
% Volga
dC2 = (V1 - 2*PriceClosed + V2) /(dV^2);
VolgaFD = 4*sqrt(V0)*(dC2*sqrt(V0) + Vega1FD/4/V0);

% Theta
T1 = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,   r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);
T2 = HestonPriceGaussLaguerre(PutCall,S0,K,Mat-dt,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);
ThetaFD = -(T1 - T2)/dt;

%% Output the results
clc;
fprintf('Stock price grid size  %5.0f\n', nS+1)
fprintf('Volatility grid size   %5.0f\n', nV+1)
fprintf('Number of time steps   %5.0f\n', nT)
fprintf('---------------------------------------------------------\n')
fprintf('Greek            Exact            PDE              FD\n')
fprintf('---------------------------------------------------------\n')
fprintf('Price         %10.5f      %10.5f            \n',PriceClosed,PricePDE)
fprintf('Delta         %10.5f      %10.5f      %10.5f\n',DeltaClosed,DeltaPDE,DeltaFD)
fprintf('Gamma         %10.5f      %10.5f      %10.5f\n',GammaClosed,GammaPDE,GammaFD)
fprintf('Vega #1       %10.5f      %10.5f      %10.5f\n',Vega1Closed,Vega1PDE,Vega1FD)
fprintf('Vanna         %10.5f      %10.5f      %10.5f\n',VannaClosed,VannaPDE,VannaFD)
fprintf('Volga         %10.5f      %10.5f      %10.5f\n',VolgaClosed,VolgaPDE,VolgaFD)
fprintf('Theta         %10.5f      %10.5f      %10.5f\n',ThetaClosed,ThetaPDE,ThetaFD)
fprintf('--------------------------------------------------------\n')

