% Heston (1993) OTM calls and puts by 32-point Gauss-Laguerre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% Uses the original Heston integrand or the Carr-Madan FFT integrand
% By Fabrice Douglas Rouah

clc; clear;

% Weights and abscissas 
[x w] = GenerateGaussLaguerre(32);

% Define the parameters and inputs
S = 1;           % Spot price.
T = 1;           % Time to maturity.
r = 0.03;        % Risk free rate.
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.25;    % Heston parameter: mean reversion level.
sigma = 0.3;     % Heston parameter: volatility of vol
v0 = .05;        % Heston parameter: initial variance.
rho = -0.8;      % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
alpha = 1.1;     % sinh damping factor

%% Price an OTM put ---------------------------------
Kp = .95;
PutCall = 'P';

% The OTM put price using Heston, Carr-Madan, damped Carr-Madan
HestonOTMPut = HestonPriceGaussLaguerre('Heston',S,Kp,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);
CarrMadanOTMPut = HestonPriceGaussLaguerre('CarrMadan',S,Kp,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);
CarrMadanOTMPutDamped = HestonPriceGaussLaguerre('CarrMadanDamped',S,Kp,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);

%% Price an OTM Call ---------------------------------
Kc = 1.05;
PutCall = 'C';

% The OTM call price using Heston, Carr-Madan, damped Carr-Madan
HestonOTMCall = HestonPriceGaussLaguerre('Heston',S,Kc,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);
CarrMadanOTMCall = HestonPriceGaussLaguerre('CarrMadan',S,Kc,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);
CarrMadanOTMCallDamped = HestonPriceGaussLaguerre('CarrMadanDamped',S,Kc,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);

%% Display the results
fprintf('\n')
fprintf('Spot = %4.2f \n',S);
fprintf('Method               OTM put K = %4.2f    OTM call K = %4.2f\n',Kp,Kc)
fprintf('----------------------------------------------------------\n')
fprintf('Heston               %10.4f           %10.4f \n', HestonOTMPut,HestonOTMCall);
fprintf('Carr Madan undamped  %10.4f           %10.4f \n', CarrMadanOTMPut,CarrMadanOTMCall);
fprintf('Carr Madan damped    %10.4f           %10.4f \n', CarrMadanOTMPutDamped,CarrMadanOTMCallDamped);
fprintf('----------------------------------------------------------\n')

clear HestonOTMPut HestonOTMCall CarrMadanOTMPut CarrMadanOTMCall
clear CarrMadanOTMPutDamped CarrMadanOTMCallDamped Kc Kp

%% Redefine the spot and strikes
S = 25;

% Price an OTM put
Kp = 20;
PutCall = 'P';
HestonOTMPut = HestonPriceGaussLaguerre('Heston',S,Kp,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);
CarrMadanOTMPut = HestonPriceGaussLaguerre('CarrMadan',S,Kp,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);

% Price an OTM Call
Kc = 30;
PutCall = 'C';
HestonOTMCall = HestonPriceGaussLaguerre('Heston',S,Kc,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);
CarrMadanOTMCall = HestonPriceGaussLaguerre('CarrMadan',S,Kc,T,r,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall,alpha);

fprintf(' \n')
fprintf('Spot = %4.2f \n',S);
fprintf('Method               OTM put K = %4.2f    OTM call K = %4.2f\n',Kp,Kc)
fprintf('----------------------------------------------------------\n')
fprintf('Heston               %10.4f           %10.4f \n', HestonOTMPut,HestonOTMCall);
fprintf('Carr Madan undamped  %10.4f           %10.4f \n', CarrMadanOTMPut,CarrMadanOTMCall);
fprintf('----------------------------------------------------------\n')
