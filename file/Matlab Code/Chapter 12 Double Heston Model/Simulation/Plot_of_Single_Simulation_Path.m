% Double Heston example from P. Gauthier and D. Possamai
% "Efficient Simulation of the Double Heston Model"

clc; clear;

% Spot price, risk free rate, dividend yield
S0 = 100;
K  = 100;
Mat = 1;
rf = 0.03;
q  = 0.02;
PutCall = 'C';

% Double Heston parameter values
v01    =  0.05;    v02    = 0.03;
sigma1 =  0.2;     sigma2 =  0.1;
rho1   = -0.9;     rho2   =  0.9;
kappa1 =  1.2;     kappa2 =  1.2;
theta1 =  0.05;    theta2 = 0.05;

% Stack the parameters into a single vector
param(1) = kappa1;
param(2) = theta1;
param(3) = sigma1;
param(4) =    v01;
param(5) =   rho1;
param(6) = kappa2;
param(7) = theta2;
param(8) = sigma2;
param(9) =    v02;
param(10)=   rho2;

% Simulation settings
NS = 1;
NT = 250;


%% Obtain true price and simulated prices
[S V1 V2 QE] = DHQuadExpSim(param,S0,K,Mat,rf,q,NT,NS,PutCall);

%% Plot the results
X = (1:NT);
[a,h1,h2] = plotyy(X,V2,X,S,'plot');
set(h1,'Color','r');
set(h2,'Color','k');
set(a(1),'YColor','k');
set(a(2),'YColor','k');
hold on
plot(X,V1,'g')
legend('Stock Price','Low Var, Pos rho','High Var, Neg rho')
hold off


