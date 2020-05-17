% Heston Local Volatility.
clc; clear;

%% Parameter estimates options with 222 days to expiry, put options on BP
param = [0.57942364171876   0.00701250929769   0.00000645051426   0.63913150446556  -0.80541286329551];
v0    = param(1);
theta = param(2); 
kappa = param(3); 
sigma = param(4); 
rho   = param(5);
lambda = 0;
T = 222 / 365;

% Risk free rate and dividend yield
rf = 0.0;
q = 0.0;

% "Little Trap" (0) or Original Heston (1) formulation
trap = 1;

%% Weights and abscissas
[x w] = GenerateGaussLaguerre(32);

% Increments for finite differences for derivatives
dt = 1e-5;
dK = 1e-1;

%% Local Volatilities
K = [15:65];
S0 = 30.67;
xT = log(K./S0);
for k=1:length(K)
	LVFD(k) = HestonLVFD(S0,K(k),T,rf,kappa,theta,sigma,lambda,v0,rho,trap,x,w,dt,dK);
	LVAN(k) = HestonLVAnalytic(S0,K(k),T,rf,kappa,theta,sigma,lambda,v0,rho,x,w,trap);
	LVAP(k) = HestonLVApprox(S0,K(k),T,kappa,theta,sigma,v0,rho);
end

%% Implied volatilities
PutCall = 'P';
a = 0.001;
b = 10;
Tol = 1e-10;
MaxIter = 1000;
for k=1:length(K)
	CallPrice = HestonCallGaussLaguerre(S0,K(k),T,rf,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
	ModelPrice(k) = CallPrice - S0 + exp(-rf*T)*K(k);
	IVm(k) = BisecBSIV(PutCall,S0,K(k),rf,T,a,b,ModelPrice(k),Tol,MaxIter);
end

%% Plot the result
plot(K,LVAP,'r-',K,LVAN,'k-',K,LVFD,'ko',K,IVm ,'rx-')
legend('Local Vol Approximate', 'Local Vol Analytic', 'Local Vol Finite Difference', 'Implied Volatility')
xlabel('Strike')
ylabel('Local and Implied Volatility')

%% Slopes
for k=2:length(K)
    slopeIV(k) = (IVm(k)  - IVm(k-1)) /(K(k-1) - K(k));
    slopeLV(k) = (LVAN(k) - LVAN(k-1))/(K(k-1) - K(k));
    ratio(k) = slopeLV(k)/slopeIV(k);
end
SlopeRatio = mean(ratio(2:length(K)));
fprintf('The ratio of Local vol to implied vol is %4.2f \n',SlopeRatio)