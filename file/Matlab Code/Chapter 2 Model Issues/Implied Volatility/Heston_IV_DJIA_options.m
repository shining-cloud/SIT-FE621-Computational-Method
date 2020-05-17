clc; clear;

%% Read in the puts IV
IV = [...
    19.62	19.47	20.19	21.15;
    19.10	19.05	19.80	20.82;
    18.60	18.61	19.43	20.57;
    18.10	18.12	19.07	20.21;
    17.61	17.74	18.71	20.00;
    17.18	17.43	18.42	19.74;
    16.71	17.06	18.13	19.50;
    16.44	16.71	17.83	19.27;
    16.45	16.41	17.60	18.99;
    16.61	16.25	17.43	18.84;
    17.01	16.02	17.26	18.62;
    17.55	16.10	17.16	18.46;
    17.96	16.57	17.24	18.42]./100;
K = (124:136);
T = [37 72 135 226]./365;
[NK NT] = size(IV);
S0 = 129.14;
rf = 0;
q  = 0;
PutCall = repmat('P',NK,NT);

%% Read in the parameter estimates and weights and abscissas

% Read in the weights and abscissas
[x w] = GenerateGaussLaguerre(32);

% Hard code the parameter estimates
param = [8.8799, 0.0674, 3.6706, 0.0435, -0.4171];
kappa = param(1);
theta = param(2);
sigma = param(3);
v0    = param(4);
rho   = param(5);
lambda = 0;

%% Fit the model prices
% "Little Trap" formulation for the Heston characteristic function
% Obtain the put prices through put call parity
trap = 1;
for k=1:NK
	for t=1:NT
		CallPrice = HestonCallGaussLaguerre(S0,K(k),T(t),rf,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
        ModelPrice(k,t) = CallPrice - S0 + exp(-rf*T(t))*K(k);
	end
end

%% Find the Model Implied Volatilities
% Starting estimates for the bisection algorithm
a = 0.001;
b = 10;
Tol = 1e-10;
MaxIter = 1000;
for k=1:NK
	for t=1:NT
		IVm(k,t) = BisecBSIV(PutCall(k,t),S0,K(k),rf,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
	end
end

%% Plot the market and model IVs
for t=1:NT
	subplot(2,2,t);
	plot(K,IV(:,t),'kx-',K,IVm(:,t),'ro-')
	legend('Market IV', 'Heston IV')
    if t==4
        legend('Market IV', 'Heston IV','Location','SouthWest')
    end
    title(['Maturity ' num2str(T(t)*365) ' days']);
    xlim([124 136]);
    ylim([.15 .22]);
	xlabel('Strike')
	ylabel('Implied Vol')
end

