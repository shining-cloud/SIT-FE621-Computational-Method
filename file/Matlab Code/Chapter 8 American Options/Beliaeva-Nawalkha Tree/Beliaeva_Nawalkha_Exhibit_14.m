% Beliaeva and Nawaklha Tree for the Heston model
% Reproduces Exhibit 14 of their paper

clc; clear;

% Heston parameters, settings, and option features
kappa = 3;
theta = 0.04;
sigma = 0.1;
rf = 0.05;
q = 0;
trap = 1;
PutCall = 'P';
Strike = 100;

% Integration points and weights for closed form price
[x w] = GenerateGaussLaguerre(32);

% Settings for the tree
NT = 15;
threshold = 0;

% Spot price, correlation, maturity, V0 parameter, and true prices
S0 = repmat([90 100 110],1,12);
rho = repmat(-[.1 .1 .1 .7 .7 .7],1,6);
T = [repmat(1/12,1,12) repmat(1/4,1,12) repmat(1/2,1,12)];
V0 = repmat([0.04.*ones(1,6) 0.16.*ones(1,6)],1,3);
BNPrice = [10 2.1254 0.1091 9.9997 2.1267 0.1274 10.71 4.2158 1.1667 10.6804 4.2140 1.1939...
    10.1706 3.4747 0.7736 10.1206 3.4807 0.8416 12.1819 6.4964 3.0914 12.1122 6.4899 3.1456 ...
    10.6478 4.6473 1.6832 10.5637 4.6636 1.7874 13.3142 8.0083 4.5454 13.2172 7.9998 4.6201];

% Loop through and find entries for Exhibit 14
fprintf(' S0    rho  sqrt(v0)   T   BNAmerPutTree  AmerPutCV    TruePrice      EuroPut\n')
fprintf('--------------------------------------------------------------------------------\n')
for k=1:36
    param = [kappa theta sigma V0(k) rho(k)]; 
    % Closed form European Put
    EuroClosed(k) = HestonPriceGaussLaguerre(PutCall,S0(k),Strike,T(k),rf,q,param,trap,x,w);
    [EuroTree(k) AmerTree(k) Euro Amer Yt V X  Prob] = BuildBivariateTree3(S0(k),PutCall,Strike,T(k),rf,NT,kappa,theta,sigma,V0(k),rho(k),threshold);
    % American put by Control Variate technique
    AmerCV(k) = AmerTree(k) + (EuroClosed(k) - EuroTree(k));
    fprintf('%4.0f %5.1f %5.1f %10.4f %10.4f %12.4f %12.4f %12.4f \n',S0(k),rho(k),sqrt(V0(k)),T(k),AmerTree(k),AmerCV(k),BNPrice(k),EuroClosed(k));
    if mod(k,12) == 0
        fprintf('--------------------------------------------------------------------------------\n')
    end
end

