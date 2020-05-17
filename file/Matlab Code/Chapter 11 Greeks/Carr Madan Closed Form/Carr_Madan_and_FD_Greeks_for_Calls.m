 % Comparison of Greeks in Carr-Madan (1999) call price
% Finite difference approximations and exact form

clc; clear;

% Gauss Laguerre abscissas and weights
[x w] = GenerateGaussLaguerre(32);

% Define the parameters and inputs
S = 100;         % Spot price.
K = 100;         % Strike
T = 0.25;        % Time to maturity.
r = 0.05;        % Risk free rate.
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.05;    % Heston parameter: mean reversion level.
sigma = 0.1;     % Heston parameter: volatility of vol
v0    = 0.05;    % Heston parameter: initial variance.
rho   = -0.9;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
alpha = 1.75;    % Carr-Madam dampening factor

%% Greeks by finite differences
% Price
PriceFD = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);

% Delta
ds = .01;
S1  = HestonCallGaussLaguerre('CarrMadan',alpha,S+ds,K,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
S0  = HestonCallGaussLaguerre('CarrMadan',alpha,S   ,K,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
S1_ = HestonCallGaussLaguerre('CarrMadan',alpha,S-ds,K,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
DeltaFD = (S1 - S1_)/(2*ds);
GammaFD = (S1 - 2*S0 + S1_)/ds^2;

% Theta
dt = 1e-4;
T1  = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T-dt,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
T1_ = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T+dt,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
ThetaFD = (T1 - T1_)/(2*dt);

% Rho
dr = 1e-5;
R1  = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T,r+dr,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
R1_ = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T,r-dr,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
RhoFD = (R1 - R1_)/(2*dr);

% Vega1
dv0 = 1e-4;
Ve1  = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T,r,kappa,theta,sigma,lambda,v0+dv0,rho,trap,x,w);
Ve1_ = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T,r,kappa,theta,sigma,lambda,v0-dv0,rho,trap,x,w);
Vega1FD = (Ve1 - Ve1_)/(2*dv0) * 2*sqrt(v0);

% Vanna
Va1  = HestonCallGaussLaguerre('CarrMadan',alpha,S+ds,K,T,r,kappa,theta,sigma,lambda,v0+dv0,rho,trap,x,w);
Va2  = HestonCallGaussLaguerre('CarrMadan',alpha,S+ds,K,T,r,kappa,theta,sigma,lambda,v0-dv0,rho,trap,x,w);
Va3  = HestonCallGaussLaguerre('CarrMadan',alpha,S-ds,K,T,r,kappa,theta,sigma,lambda,v0+dv0,rho,trap,x,w);
Va4  = HestonCallGaussLaguerre('CarrMadan',alpha,S-ds,K,T,r,kappa,theta,sigma,lambda,v0-dv0,rho,trap,x,w);
VannaFD = (Va1 - Va2 - Va3 + Va4)/4/dv0/ds*2*sqrt(v0);

% Volga
Ve0  = HestonCallGaussLaguerre('CarrMadan',alpha,S,K,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
dC2 = (Ve1 - 2*Ve0 + Ve1_)/dv0^2;
VolgaFD = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1FD/4/v0);


%% Greeks by closed form
[x w] = GenerateGaussLaguerre(32);
Price = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Price',x,w);
Delta = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Delta',x,w);
Gamma = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Gamma',x,w);
Theta = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Theta',x,w);
Rho   = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Rho',x,w);
Vega1 = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Vega1',x,w);
Vanna = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Vanna',x,w);
Volga = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,'Volga',x,w);


%% Write the results
fprintf('Carr-Madan Greeks by finite differences and by closed form\n');
fprintf('-------------------------------------\n')
fprintf('Greek            F.D.       Closed \n')
fprintf('-------------------------------------\n')
fprintf('Call Price   %10.4f  %10.4f \n',PriceFD,Price);
fprintf('Delta        %10.4f  %10.4f \n',DeltaFD,Delta);
fprintf('Gamma        %10.4f  %10.4f \n',GammaFD,Gamma);
fprintf('Theta        %10.4f  %10.4f \n',ThetaFD,Theta);
fprintf('Rho          %10.4f  %10.4f \n',RhoFD,Rho);
fprintf('Vega1        %10.4f  %10.4f \n',Vega1FD,Vega1);
fprintf('Vanna        %10.4f  %10.4f \n',VannaFD,Vanna);
fprintf('Volga        %10.4f  %10.4f \n',VolgaFD,Volga);
fprintf('-------------------------------------\n')

