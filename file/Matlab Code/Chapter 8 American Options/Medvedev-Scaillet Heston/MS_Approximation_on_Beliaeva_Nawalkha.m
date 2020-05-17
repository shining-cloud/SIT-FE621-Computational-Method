% Medvedev-Scaillet Approximation to American put prices 
% of Beliaeva and Nawalkha

clc; clear;

% Option price settings
Strike = 100;
T = 0.5;
r = 0.05;
q = 0;

% Heston parameters
kappaV = 3;
thetaV = 0.04;
sigmaV = 0.1;
lambda = 0; 
trap = 1;

% Spot prices
S = repmat([90 100 110],1,12)';
% Correlations
rho = repmat([-.1 -.1 -.1 -.7 -.7 -.7],1,6)';
% Variances
v0 = repmat([(0.2^2).*ones(1,6) (0.4^2).*ones(1,6)],1,3)';
% Maturirites
T = [(1/12).*ones(1,12) (0.25).*ones(1,12) (0.5).*ones(1,12)]';
% True prices
True = [10.0000 2.1254 0.1091  9.9997 2.1267 0.1274 10.7100 4.2158 1.1667 10.6804 4.2140 1.1939 ...
        10.1706 3.4747 0.7736 10.1206 3.4807 0.8416 12.1819 6.4964 3.0914 12.1122 6.4899 3.1456 ...
        10.6478 4.6473 1.6832 10.5637 4.6636 1.7874 13.3142 8.0083 4.5454 13.2172 7.9998 4.6201]';

% Settings for Simpson's rule
A = 1e-20;
B = 100;
N = 1000;
method = 3;

% Infinite barrier for the European put
yinf = 10000;

% Select the number of terms in the expansion (3,4 or 5)
NumTerms = 4;
%MSPut = @MSPutHeston;


%% Find the Medvedev-Scaillet 5-term approximation
for k=1:length(S);
    params = [kappaV thetaV sigmaV v0(k) rho(k)];
    [EuroPutClosed AmerPutMS AmerPut(k) EEP(k) theta(k) y(k)] = MSPrice(S(k),Strike,T(k),r,q,params,trap,method,A,B,N,NumTerms,yinf);
end

error = True - AmerPut';
TotalError = sum(abs(error));

%% Display the results
clc;
fprintf(['Medvedev-Scaillet ' num2str(NumTerms) '-term approximation \n'])
fprintf('Beliaeva and Nawalkha prices \n')
fprintf('----------------------------------------------------------- \n')
fprintf(' S(0)  TruePrice   AmerPut    $Error      Barrier     EEP\n')
fprintf('----------------------------------------------------------- \n')
for k=1:length(S)
	fprintf('%4.0f %10.4f %10.4f %10.4f %10.4f %10.4f\n',S(k),True(k),AmerPut(k),error(k),y(k),EEP(k));
end
fprintf('----------------------------------------------------------- \n')
fprintf('Total Absolute Error %6.4f \n',TotalError)
fprintf('----------------------------------------------------------- \n')
