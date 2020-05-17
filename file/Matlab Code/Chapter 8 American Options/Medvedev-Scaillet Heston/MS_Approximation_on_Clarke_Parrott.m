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
N = 2000;
method = 3;

% Infinite barrier for the European put
yinf = 1e4;

% Select the number of terms in the expansion (3,4 or 5)
NumTerms = 4;

%% Find the Medvedev-Scaillet Heston price
for k=1:5
    [EuroPutClosed(k) AmerPutMS AmerPut(k) EEP theta(k) y(k)] = MSPrice(S(k),Strike,T,r,q,params,trap,method,A,B,N,NumTerms,yinf);
end

% Errors with 5 terms
error = TruePrice - AmerPut;
TotalError = sum(abs(error));


%% Write the results
fprintf(['Medvedev-Scaillet ' num2str(NumTerms) '-term approximation \n'])
fprintf('Clarke and Parrott prices \n')
fprintf('-----------------------------------------------------------\n');
fprintf('Spot  TruePrice    Approx     Error     Barrier   Moneyness\n');
fprintf('-----------------------------------------------------------\n');
for k=1:5
	fprintf('%3.0f %10.4f %10.4f %10.4f %10.4f %10.4f\n',...
             S(k),TruePrice(k),AmerPut(k),error(k),y(k),theta(k));
end
fprintf('-----------------------------------------------------------\n');
fprintf(['Total Absolute Error with ' num2str(NumTerms) ' terms %10.8f \n'],TotalError)
fprintf('-----------------------------------------------------------\n');

 