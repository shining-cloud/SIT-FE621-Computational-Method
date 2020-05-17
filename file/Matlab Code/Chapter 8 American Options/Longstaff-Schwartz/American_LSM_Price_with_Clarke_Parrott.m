% Evaluation of American put options fromt the Heston model.
% Benchmark against values in Clarke and Parrot (1999)

clc; clear;

% Import the spot prices and true values
S = [8 9 10 11 12];
TruePrice = [2.0000, 1.107641, 0.520030, 0.213668, 0.082036];

% Settings from Clarke and Parrot (1999)
K = 10;
kappa = 5;
theta = 0.16;
sigma = 0.9;
v0 = 0.0625;
rho = 0.1;
lambda = 0;
T = 1/4;
r = 0.1;
q = 0.0;
params = [kappa theta sigma v0 rho lambda];

% Settings for the European prices
PutCall = 'P';
[x w] = GenerateGaussLaguerre(10);
trap = 1;

% Settings for the LSM algorithm
% Number of time steps and number of stock price paths
NT = 100;
NS = 1000;

% Design matrix for the LSM algorithm
XmatrixHandle = {@(y)ones(length(y),1), @(y)(1-y),@(y)1./2.*(2-4.*y-y.^2)};

% Matrices for the correlated N(0,1) random variables
Zv = randn(NT,NS);
Zs = rho.*Zv + sqrt(1-rho.^2).*randn(NT,NS);

% Loop through each stock price, calculating the prices
for s=1:length(S)
    % Generate the paths and LSM prices
    [Spaths V] = MMSim(params,S(s),T,r,q,NT,NS,Zv,Zs);
    [LSMEuro(s) LSMAmer(s)] = LSM(Spaths',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    % Generate the European price and errors
    ClosedEuro(s) = HestonPriceGaussLaguerre(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
    % Control variate American price
    CV(s) = ClosedEuro(s) + (LSMAmer(s) - LSMEuro(s));
end

% The errors 
Error   = TruePrice - LSMAmer;
ErrorCV = TruePrice - CV;
ErrorE  = ClosedEuro - LSMEuro;
TotalError   = sum(abs(Error));
TotalErrorCV = sum(abs(ErrorCV));
TotalErrorE  = sum(abs(ErrorE));

%% Display the results for American puts
disp('Comparison of Clarke and Parrott American Put prices with Least-Squares Monte-Carlo')
disp(['LSM uses ' num2str(NT) ' time steps, and ' num2str(NS) ' stock paths '])
disp('-----------------------------------------------------------------------------')
disp('  S(0)   TruePrice   LSMAmer  $Error   ControlVar  $Error   TrueEuro  LSMEuro')
disp('-----------------------------------------------------------------------------')
disp(num2str([S' TruePrice' LSMAmer' Error' CV' ErrorCV' ClosedEuro' LSMEuro'],'%10.4f'))
disp('-----------------------------------------------------------------------------')
fprintf('AbsError %28.4f %19.4f %19.4f \n',TotalError,TotalErrorCV,TotalErrorE')
disp('-----------------------------------------------------------------------------')


