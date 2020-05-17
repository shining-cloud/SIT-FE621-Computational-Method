% Illustration of Kahl-Jackel (2005) integration range

clc; clear;

S  = 10;         % Stock price
K  = 7;          % Strike price
rf = 0.06;       % Risk free rate
q  = 0.04;       % Dividend yield
kappa = 1;       % Mean reversion speed
theta = 0.06;    % Mean reversion level
sigma = 0.5;     % Volatility of variance
v0  = 0.06;      % Spot (initial) variance
rho = -0.8;      % Correlation
T = 1/12;        % Maturity in years
PutCall = 'C';   % 'P'ut or 'C'all
trap = 1;        % "Little Trap" formulation
lambda = 0;      % Risk parameter


%% Compare the call values under different integration schemes
% 32-point Gauss Laguerre integration, Heston integrand
param = [kappa theta sigma v0 rho];
[x w] = GenerateGaussLaguerre(32);
PriceGLa = HestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);

% 32-point Gauss-Lobatto integration, Heston integrand
[x w] = GenerateGaussLobatto(32);
A = 1e-10;
B = 100;
PriceGLo = HestonPriceGaussLegendre(PutCall,S,K,T,rf,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,A,B);

% 32-point Gauss-Lobatto Integration, Kahl Jackel integrand
PriceKJ = HestonPriceKahlJackel(PutCall,S,K,T,rf,q,param,trap,x,w);

%% Display the results
disp(['32 point Gauss Laguerre, Heston integrand       ' num2str(PriceGLa)])
disp(['32 point Gauss Lobatto, Heston integrand        ' num2str(PriceGLo)])
disp(['32 point Gauss Lobatto, Kahl-Jackel integrand   ' num2str(PriceKJ)])


%% Analysis of the Kahl-Jackel Integrand
N = 20000;       % Number of points in integrand range
du = 1/(N-1);    % Increment of integrand range
U = [0:du:1];

% Integrand at continuous values
for k=1:N;
	y(k) = HestonIntegrandKahlJackel(U(k),PutCall,S,K,T,rf,q,param,trap);
end

% Integrand at Gauss-Lobatto abscissas only
for k=1:length(x)
	phi = 1/2*x(k) + 1/2;
	z(k) = HestonIntegrandKahlJackel(phi,PutCall,S,K,T,rf,q,param,trap);
end

%% Plot the integrand
plot(U,y,'k-',0.5.*x +1/2,z,'ro')
legend('Kahl-Jackel Integrand', 'K-J Integrand at Gauss-Lobatto abscissas','Location','NorthWest')


