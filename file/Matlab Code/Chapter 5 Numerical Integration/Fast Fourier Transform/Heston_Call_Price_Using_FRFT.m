% Fractional Fast Fourier Transform for the Heston Model
 
clc; clear;

% Required inputs
N = 2^12;          % Number points for the FFT
S0 = 50;           % Spot price.
r = 0.03;          % Risk free rate.
q = 0.05;           % Dividend yield
tau = .5;          % Time to maturity.
kappa  = 0.2;      % Heston parameter: mean reversion speed.
theta  = 0.05;     % Heston parameter: mean reversion level.
sigma  = 0.3;      % Heston parameter: volatility of vol
lambda = 0;        % Heston parameter: risk.
rho    = -0.7;     % Heston parameter: correlation
v0     = 0.05;     % Heston parameter: initial variance.
trap = 1;          % Heston trap (1) or original Heston (0) formulation
alpha = 1.5;       % Dampening factor
uplimit = 100;     % Upper limit of integration
fast = 1;          % Choice of fast (1) or slow (0) algorithm

% 32-point Gauss Laguerre abscissas and weights
[x w] = GenerateGaussLaguerre(32);

%% Obtain the FRFT call prices, trapezoidal rule
% Integration increment
eta = 1e-2;
lambdainc = eta;

% Alternative: Select lambdainc so that the strike range is (Spot +/ 5);
lambdainc = 2/N*log(S0/(S0-30));

rule = 'T';
[CallT K lambdainc eta] = HestonCallFRFT(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,trap,alpha,rule,eta,lambdainc);

% Verify the strike range
disp(['FRFT strike range from ' num2str(K(1)) ' to ' num2str(K(end))])
disp(' ');

% Find the ATM prices and obtain the exact prices
%ATM = [find(round(K)==S0)-2:find(round(K)==S0)+4];
ATM = [find(round(K)==S0)+19:find(round(K)==S0)+25];     % Alternative
CallT = CallT(ATM);


%% Obtain the FRFT call prices, Simpson's rule
% Integration increment and strike increment
rule = 'S';
[CallS K lambdainc eta] = HestonCallFRFT(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,trap,alpha,rule,eta,lambdainc);

% Find the ATM prices and obtain the exact prices
CallS = CallS(ATM);

%% Fit the exact prices and calculate the error
K = K(ATM);
for j=1:length(K);
	CallE(j) = HestonPriceGaussLaguerre('C',S0,K(j),tau,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
end
ErrorT = (CallT-CallE')./CallE'.*100;
ErrorS = (CallS-CallE')./CallE'.*100;

%% Print the results
disp('Strike       Exact       FRFT Trapz  FRFT Simp  %Error Tr   %Error Sim')
disp('----------------------------------------------------------------------')
disp(num2str([K CallE' CallT CallS ErrorT ErrorS],'%12.4f'))
disp(' ')

% Print the increments for integration and for log strikes
disp(['Integration increment ' num2str(eta,'%10.6f')])
disp(['Log strike increment  ' num2str(lambdainc,'%16.6f\n')])


