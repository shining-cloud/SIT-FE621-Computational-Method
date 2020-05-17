% Heston (1993) Call price by Newton-Coates Formulas
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% By Fabrice Douglas Rouah

clc; clear;

% Define the parameters and inputs
S = 100;         % Spot price.
T = 0.5;         % Time to maturity.
r = 0.03;        % Risk free rate.
q = 0.02;        % Dividend yield
kappa = 0.2;     % Heston parameter: mean reversion speed.
theta = 0.25;    % Heston parameter: mean reversion level.
sigma = 0.3;     % Heston parameter: volatility of vol
lambda = 0;      % Heston parameter: risk.
v0 = .02;        % Heston parameter: initial variance.
rho = -0.8;      % Heston parameter: correlation
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
method = 1;      % 1 = mid-point rule
                 % 2 = trapezoidal rule
				 % 3 = Simpson's rule
				 % 4 = Simpson's 3/8 rule
a = 1e-8;        % Lower limit for Newton-Cotes
b = 100;         % Upper limit for Newton-Cotes
N = 5000;        % Number of abscissas for Newton-Cotes

Method = {'Mid-Point   ' 
	      'Trapezoidal ' 
		  'Simpsons    ' 
		  'Simpsons 3/8'};

K = [80:10:120];

for k=1:length(K);
	disp(['------ ' 'Strike = ' num2str(K(k)) ' ----------------------']);
	for method=1:4;
		Call = HestonPriceNewtonCoates('C',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);
		Put  = HestonPriceNewtonCoates('P',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);
		disp([char(Method(method)) ' Call = ' num2str(Call,'%10.6f') '  Put = ' num2str(Put,'%10.6f') ])
	end
end

