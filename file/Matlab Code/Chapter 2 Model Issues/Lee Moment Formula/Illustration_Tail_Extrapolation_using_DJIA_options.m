% Extrapolation of Heston tails
% Uses DIA puts on May 10, 2012

clc; clear;

% Function for the Black Scholes put
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));
  
%% Load the DIA data
PutIV = [19.62, 19.10, 18.60, 18.10, 17.61, 17.18, 16.71, ...
         16.44, 16.45, 16.61, 17.01, 17.55, 17.96]'./100;
K = (124:136);
T = 37/365;
S = 129.14;
F = S;
M = log(K./F);
rf = 0.0;
q  = 0.0;
trap = 1;

% Parameter estimates for the single maturity : kappa theta sigma v0 rho
kappa =  4.5150;
theta =  0.1250;
sigma =  2.1149;
v0    =  0.0200;
rho   = -0.2741;
param = [kappa theta sigma v0 rho];

% Gauss Laguerre Weights and abscissas
[x w] = GenerateGaussLaguerre(32);
trap = 1;

%% Fit the model implied volatilities from the model prices
% Settings for the bisection algorith
a = 0.01;
b = 4.0;
MaxIter = 2000;
Tol = 1e-6;

% Extrapolate +/- 6 on the strikes and fit Heston prices, extract IV
dk = 4;
dK = 10;
K1 = [K(1)-dk:K(end)+dk];
for k=1:length(K1)
    SinglePrice(k) = HestonPriceGaussLaguerre('P',S,K1(k),T,rf,q,param,trap,x,w);
    SingleIV(k) = BisecBSIV('P',S,K1(k),rf,q,T,a,b,SinglePrice(k),Tol,MaxIter);
end
F = S;
M1 = log(K1./F);

%% Display the parameter estimates
fprintf('Heston model parameter estimates\n')
fprintf('kappa     theta     sigma      v0       rho\n');
fprintf('----------------------------------------------\n');
disp([num2str(param,'%10.4f')])
fprintf('----------------------------------------------\n');
fprintf(' \n')

%% Roger Lee slope bounds
HiLimit =  10;
LoLimit = -10;
[bR bL LowerAP UpperAP LowerCF UpperCF] = FindLeeBounds(S,rf,q,T,param,trap,LoLimit,HiLimit);
ptilde =  UpperAP - 1;
qtilde = -LowerAP;
fprintf('Moment bounds            Lower      Upper\n');
fprintf('------------------------------------------\n');
fprintf('Andersen & Piterbarg %10.4f %10.4f\n',LowerAP,UpperAP);
fprintf('Roger Lee q~ and p~  %10.4f %10.4f\n',qtilde,ptilde);
fprintf(' \n');

fprintf('Roger Lee slope limits\n');
fprintf('Left Slope   Right Slope \n');
fprintf('-------------------------\n');
fprintf('%8.4f  %10.4f \n',bL,bR);
fprintf('-------------------------\n');


%% Right tail extrapolation +5 strikes
K2 = [K1(end):K1(end)+dK];
M2 = [log(K2(1)/F) : .01 : log(K2(end)/F)];
increment = M2(2) - M2(1);
RightTail(1) = SingleIV(end)^2;

% Right slope
for j=2:length(M2);
    RightTail(j) = RightTail(j-1) + bR*increment;
end
RightTail = sqrt(RightTail);

%% Left tail extrapolation -5 strikes
K3 = [K1(1)-dK:K1(1)];
M3 = fliplr([log(K3(end)/F) : -0.01 : log(K3(1)/F)]);
increment = M3(2) - M3(1);
LeftTail = zeros(1,length(M3));
LeftTail(end) = SingleIV(1)^2;

% Left slope
for j=length(M3)-1:-1:1;
    LeftTail(j) = LeftTail(j+1) + bL*increment;
end
LeftTail = sqrt(LeftTail);

%% Strike ranges
fprintf(' \n');
fprintf('Strike range expansion \n');
fprintf('Range              Strikes           Moneyness \n');
fprintf('-----------------------------------------------------\n');
fprintf('In Sample     %5.0f %10.0f %10.4f %10.4f \n',K(1),K(end),M(1),M(end))
fprintf('Out-of-sample %5.0f %10.0f %10.4f %10.4f \n',K1(1),K1(end),M1(1),M1(end))
fprintf('Left tail     %5.0f %10.0f %10.4f %10.4f \n',K3(1),K3(end),M3(1),M3(end))
fprintf('Right tail    %5.0f %10.0f %10.4f %10.4f \n',K2(1),K2(end),M2(1),M2(end))
fprintf('-----------------------------------------------------\n');

%% Plot the final result
plot(M,PutIV,'ko',M1,SingleIV,'rx-',M2,RightTail,'k--',M3,LeftTail,'k--')
legend('Market IV', 'Heston IV', 'Extrapolated Tail', 'Extrapolated Tail')
xlabel('Log Moneyness')
ylabel('Implied Volatility')
axis([-0.2 0.2 0.16 0.28]) 
