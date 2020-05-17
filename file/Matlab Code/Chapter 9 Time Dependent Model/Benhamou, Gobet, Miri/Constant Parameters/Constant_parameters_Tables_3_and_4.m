% Reproduces Tables 3 and 4 of Benhamou, Gobet, and Miri
% "Time Dependent Heston Model"
% SIAM Journal on Financial Mathematics, Vol. 1, (2010)

clc; clear;

% Spot price, risk free rate, dividend yield
% Note the typo in BGM: prices in Tables 3 and 4 are from calls, not puts
S  = 100;
rf = 0.0;
q  = 0.0;
PutCall = 'C';

% Heston model parameters.
% BGM Tables 3 and 4
kappa  = 3;       % Volatility reversion speed
theta  = 0.06;    % Volatility reversion level
sigma  = 0.3;     % Volatility of variance
rho    = 0;       % Correlation
v0     = 0.04;    % Initial variance
lambda = 0;       % Risk
trap   = 1;       % Little trap formulation

% Vectorize the Heston parameters
params(1) = kappa;
params(2) = theta;
params(3) = sigma;
params(4) = rho;
params(5) = v0;
params(6) = lambda;

% Input the maturities and the strike prices
T = [.25 .5 1 2 3 5 7 10];
T = repmat(T,8,1)';
K = [70 80 90 100 110 120 125 130;
	 60 70 80 100 110 130 140 150;
	 50 60 80 100 120 150 170 180;
	 40 50 70 100 130 180 210 240;
	 30 40 60 100 140 200 250 290;
	 20 30 60 100 150 250 320 400;
	 10 30 50 100 170 300 410 520;
	 10 20 50 100 190 370 550 730];
[x w] = GenerateGaussLaguerre(32);

%% Calculate the call prices for Table 4
for j=1:8
	for i=1:8
		ApproxCall(i,j) = BGMApproxPrice(params,S,K(i,j),rf,q,T(i,j),trap,PutCall);
		ExactCall(i,j)  = HestonPriceGaussLaguerre(PutCall,S,K(i,j),T(i,j),rf,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
	end
end

fprintf('Table 4 of BGM (2010)\n')
fprintf('--------------------------------------------------------------------------------\n');
fmtstring = [repmat('%10.4f',1,7) '%10.4f\n'];
for j=1:8;
    fprintf(fmtstring,ExactCall(j,:));
    fprintf(fmtstring,ApproxCall(j,:));
end
fprintf('\n');

%% Calculate the implied vols for Table 3
a = 1e-2;
b = 5;
MaxIter = 3000;
Tol = 1e-8;
for j=1:8
	for i=1:8
		ApproxIV(i,j) = BisecBSIV(PutCall,S,K(i,j),rf,q,T(i,j),a,b,ApproxCall(i,j),Tol,MaxIter);
		ExactIV(i,j)  = BisecBSIV(PutCall,S,K(i,j),rf,q,T(i,j),a,b,ExactCall(i,j) ,Tol,MaxIter);
	end
end

fprintf('Table 3 of BGM (2010)\n')
fprintf('--------------------------------------------------------------------------------\n');
fmtstring = [repmat('%10.2f',1,7) '%10.2f\n'];
for j=1:8;
    fprintf(fmtstring,ExactIV(j,:).*100);
    fprintf(fmtstring,ApproxIV(j,:).*100);
end
fprintf('\n');

%% Output only the ATM values
disp('Table 4 - ATM values only')
fprintf('-----------------------------------\n');
disp(num2str([ExactIV(:,4).*100 ApproxIV(:,4).*100 ExactCall(:,4) ApproxCall(:,4)],'%10.2f'))
