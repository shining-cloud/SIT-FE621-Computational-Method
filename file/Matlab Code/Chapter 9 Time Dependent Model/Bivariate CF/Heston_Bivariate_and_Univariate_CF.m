% Demonstrate that prices are identical when the univariate characteristic
% function is used, or when the bivariate characteristic function is used

clc; clear;

T = 0.25;         % Time to maturity.
S = 100;          % Spot price.
r = 0.03;         % Risk free rate.
q = 0.07;          % Dividend yield
kappa  =  0.2;    % Heston parameter: mean reversion speed.
theta  =  0.05;   % Heston parameter: mean reversion level.
sigma  =  0.3;    % Heston parameter: volatility of vol
rho    = -0.8;    % Heston parameter: correlation
lambda =  0;      % Heston parameter: risk
v0     =  0.22;   % Heston parameter: initial variance.
PutCall = 'C';    % 'C'all or 'P'ut
trap = 1;         % 1 = "Little Trap" formulation
                  % 0 = Original Heston formulation
				  
[x w] = GenerateGaussLaguerre(32);
K = [90:110];
for k=1:length(K)
	% Prices using the original Heston formulation
	trap=0;
	% Univariate c.f.
	Price1(k) = HestonPriceGaussLaguerre(PutCall,S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,1);
	% Bivariate c.f.
	Price2(k) = HestonPriceGaussLaguerre(PutCall,S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,2);
	% Prices using the "Little Trap" formulation
	trap=1;
	% Univariate c.f.
	Price3(k) = HestonPriceGaussLaguerre(PutCall,S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,1);
	% Bivariate c.f.
	Price4(k) = HestonPriceGaussLaguerre(PutCall,S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,2);
end

fprintf('Equivalence of univariate and bivariate CF on prices \n')
fprintf('----------------------------------------------------\n')
fprintf('          Original Heston CF  |   Little Trap CF\n')
fprintf(' Strike  Univariate Bivariate | Univariate Bivariate\n')
fprintf('----------------------------------------------------\n')
for k=1:length(K)
    fprintf(' %5.0f %10.4f %10.4f %10.4f %10.4f\n',K(k),Price1(k),Price2(k),Price3(k),Price4(k));
end
fprintf('----------------------------------------------------\n')
