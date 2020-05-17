% Fast Fourier Transform for the Heston Model

% --------------------------------------
% INPUTS
%   N  = number of discretization points
%   uplimit = Upper limit of integration
%   S0 = spot price
%   r = risk free rate
%   tau = maturity
%   sigma = volatility
%   alpha = dampening factor
%   fast = fast versus slow algorithm.
%     fast = 1 fast version.  Uses vectorization.
%     fast = 0 slow version.  Uses loops
% --------------------------------------
% OUTPUTS
%   CallFFT = Black Scholes call prices using FFT
%   CallBS  = Black Scholes call prices using closed form
%         K = Strike prices
%       eta = increment for integration range
%    lambda = increment for log-strike range
% --------------------------------------
 
clc; clear;

% Required inputs
N = 2^10;          % Number points for the FFT
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

%% Run the Fast Fourier Transform using the trapezoidal rule
[CallT K lambdainc eta] = HestonCallFFT(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,trap,alpha,fast,'T');

% Run the Fast Fourier Transform using Simpson's rule
[CallS K lambdainc eta] = HestonCallFFT(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,trap,alpha,fast,'S');

% Obtain the results near the ATM strikes
ATM = [find(round(K)==S0)-3:find(round(K)==S0)+3];

% Truncate the outputted calls and strikes
CallT = CallT(ATM);
CallS = CallS(ATM);
K     = K(ATM);

%% Closed form price and errors
[x w] = GenerateGaussLaguerre(32);
for k=1:length(K);
	CallE(k) = HestonPriceGaussLaguerre('C',S0,K(k),tau,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
end
ErrorT = (CallT-CallE')./CallE'.*100;
ErrorS = (CallS-CallE')./CallE'.*100;

MSET = mean(abs(ErrorT));
MSES = mean(abs(ErrorS));

%% Print the results
disp('Strike       Exact       FFT Trapz   FFT Simp   %Error Tr   %Error Sim')
disp('----------------------------------------------------------------------')
disp(num2str([K CallE' CallT CallS ErrorT ErrorS],'%12.4f'))
disp(' ')

% Print the increments for integration and for log strikes
disp(['Integration increment ' num2str(eta,'%10.6f')])
disp(['Log strike increment  ' num2str(lambdainc,'%10.6f\n')])

% Print the errors
disp(['Trapezoidal FFT mean absolute error ' num2str(MSET)])
disp(['Simpsons FFT mean absolute error    ' num2str(MSES)])

