% Heston American put prices of Clarke and Parrott (1999)
% Uses Medvedev-Scaillet Expansion

clc; clear;

% Settings from Clarke and Parrott
S = [8 9 10 11 12];
TruePrice = [2.00 1.107641 0.520030 0.213668 0.082036];
Strike = 10;
r = 0.1;
q = 0;
T = 0.25;

% Heston parameters
kappaV = 5;
thetaV = 0.16;
sigmaV = 0.9;
rho = 0.1;
v0  = 0.0625;
lambda = 0;
trap = 1;
params = [kappaV thetaV sigmaV v0 rho];

% Settings for Simpson's Rule
A = 1e-20;
B = 100;
N = 10000;
method = 3;

% Settings for the optimzation
hi = 3;       % Upper limit for fminbnd
yinf = 1e4;   % Infinite barrier for European put

% Select the number of terms in the expansion (3,4 or 5)
NumTerms = 4;

%% Find the finite difference prices and Greeks
for k=1:5
    Price(k) = MSGreeksFD(params,S(k),Strike,r,q,T,method,A,B,N,yinf,hi,NumTerms,'price');
    Delta(k) = MSGreeksFD(params,S(k),Strike,r,q,T,method,A,B,N,yinf,hi,NumTerms,'delta');
    Gamma(k) = MSGreeksFD(params,S(k),Strike,r,q,T,method,A,B,N,yinf,hi,NumTerms,'gamma');
    Vega1(k) = MSGreeksFD(params,S(k),Strike,r,q,T,method,A,B,N,yinf,hi,NumTerms,'vega1');
    Vanna(k) = MSGreeksFD(params,S(k),Strike,r,q,T,method,A,B,N,yinf,hi,NumTerms,'vanna');
    Volga(k) = MSGreeksFD(params,S(k),Strike,r,q,T,method,A,B,N,yinf,hi,NumTerms,'volga');
    Theta(k) = MSGreeksFD(params,S(k),Strike,r,q,T,method,A,B,N,yinf,hi,NumTerms,'theta');
end

%% Write the results
fprintf(['Greeks using Medvedev-Scaillet ' num2str(NumTerms) '-term approximation \n'])
fprintf('Clarke and Parrott prices \n')
fprintf('-------------------------------------------------------------------------------------------\n');
fprintf('Spot  TruePrice    Approx      Delta      Gamma     Vega1      Vanna      Volga      Theta\n');
fprintf('-------------------------------------------------------------------------------------------\n');
for k=1:5
	fprintf('%3.0f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',...
             S(k),TruePrice(k),Price(k),Delta(k),Gamma(k),Vega1(k),Vanna(k),Volga(k),Theta(k));
end
fprintf('-------------------------------------------------------------------------------------------\n');



