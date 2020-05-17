% Solving the PDE for the Heston model for a European Call.
% In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
% in the Heston Modelo with Correlation" Int J of Num Analysis and
% Modeling, 2010.

clc; clear;

% Strike price, risk free rate, dividend yield, and maturity
S0  = 101.52;
K   = 100;
r   = 0.02;
q   = 0.05;
Mat = 0.15;

% Heston parameters.  Case 1 of Hout and Foulon (Table 1)
kappa =  1.5;
theta =  0.04;
sigma =  0.3;
rho   = -0.9;
V0    =  0.05412;
lambda = 0;

% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0;   Smax = 2*K;
Vmin = 0;   Vmax = 0.5;
Tmin = 0;   Tmax = Mat;

% Number of grid points for U(S,v,t) is nS+1, nV+1, nT+1
nS = 40;        % Stock price
nV = 40;        % Volatility
nT = 999;      % Maturity

% Factor used to round the option values U
RoundFactor = 100000;

%% Decreasing error with increasing Stock price and volatility grid sizes
params = [kappa theta sigma V0 rho lambda];
trap = 1;
PutCall = 'C';
numS = 19:5:nS-1;
numV =  9:10:nV;

% Closed form price
[x w] = GenerateGaussLaguerre(32);
TruePrice = HestonPriceGaussLaguerre(PutCall,S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);

%% Uniform Grid -- Loop through and find the PDE prices and the error
for s=1:length(numS)
	for v=1:length(numV)
		% Build the grids
		nS = numS(s);
		nV = numV(v);
		% Increment for Stock Price, Volatility, and Maturity
		ds = (Smax-Smin)/nS;
		dv = (Vmax-Vmin)/nV;
		dt = (Tmax-Tmin)/nT;
		% Number of points on the 3-D grid
		NS = nS+1;        % Stock price
		NV = nV+1;        % Volatility
		NT = nT+1;        % Maturity
		% Vectors for the Stock Price, Volatility, and Maturity
		S = [0:NS-1].*ds;
		V = [0:NV-1].*dv;
		T = [0:NT-1].*dt;
		U = HestonExplicitPDE(params,K,r,q,S,V,T);
		PriceU(s,v) = interp2(V,S,U,V0,S0);
		ErrorU(s,v) = PriceU(s,v) - TruePrice;
        fprintf('Completed uniform grid price for NS = %3.0f and NV = %3.0f \n',NS,NV);
	end
end


%% Non-Uniform Grid -- Loop through and find the PDE prices and the error
clear S V

% The maturity time increment and grid
dt = (Tmax-Tmin)/nT;
T = [0:NT-1].*dt;

for s=1:length(numS)
	nS = numS(s);
	for v=1:length(numV)
		nV = numV(v);
		% Number of points on the 3-D grid
		NS = nS+1;
		NV = nV+1;
		NT = nT+1;
		% The stock price grid S(i)
		c = K/5;        % Value used by Int 'T Hout and Foulon
		dz = 1/(NS-1)*(asinh((Smax-K)/c) - asinh(-K/c));
		for i=1:NS;
			z(i) = asinh(-K/c) + (i-1)*dz;
			S(i) = K + c*sinh(z(i));
		end
		% The volatility grid V(j)
		d = Vmax/500;   % Value used by Int 'T Hout and Foulon
		dn = asinh(Vmax/d)/(NV-1);
		for j=1:NV
			n(j) = (j-1)*dn;
			V(j) = d*sinh(n(j));
		end
		% Solve the PDE using the (s,v)th grid
		U = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T);
		PriceN(s,v) = interp2(V,S,U,V0,S0);
		ErrorN(s,v) = PriceN(s,v) - TruePrice;
        fprintf('Completed non-uniform grid price for NS = %4.0f and NV = %4.0f \n',NS,NV);
		clear V
	end
	clear S
end


%% Output the results
clc;
fprintf('-------------------------------------------------------------------\n');
fprintf('                               Uniform Grid       Non-Uniform Grid \n');
fprintf(' NS  NV   NT    ExactPrice   PDEPrice   Error     PDEPrice    Error\n');
fprintf('-------------------------------------------------------------------\n');
for s=1:length(numS)
	for v=1:length(numV)
        fprintf('%3.0f %3.0f %5.0f %10.4f %10.4f %10.4f %10.4f %10.4f\n', numS(s)+1,numV(v)+1,NT,TruePrice,PriceU(s,v),ErrorU(s,v),PriceN(s,v),ErrorN(s,v));
	end
end
fprintf('-------------------------------------------------------------------\n');


%% Plot the prices
True = TruePrice.*ones(1,5);
plot(numS+1,PriceU(:,1),'ko:',numS+1,PriceU(:,2),'ko--',numS+1,PriceU(:,3),'ko-.',numS+1,PriceU(:,4),'ko-',...
     numS+1,PriceN(:,1),'ro:',numS+1,PriceN(:,2),'ro--',numS+1,PriceN(:,3),'ro-.',numS+1,PriceN(:,4),'ro-',numS+1,True);
legend('Uniform NV=10',    'Uniform NV=20',    'Uniform NV=30',    'Uniform NV=40',...
       'Non-Uniform NV=10','Non-Uniform NV=20','Non-Uniform NV=30','Non-Uniform NV=40','True Price');
xlabel('Number of stock price points');


% %% Plot the errors
% plot(numS,ErrorU(:,1),'k:',numS,ErrorU(:,2),'k--',numS,ErrorU(:,3),'k-.',numS,ErrorU(:,4),'k-',...
%      numS,ErrorN(:,1),'r:',numS,ErrorN(:,2),'r--',numS,ErrorN(:,3),'r-.',numS,ErrorN(:,4),'r-');
% legend('Uniform NV=10',    'Uniform NV=20',    'Uniform NV=30',    'Uniform NV=40',...
%        'Non-Uniform NV=10','Non-Uniform NV=20','Non-Uniform NV=30','Non-Uniform NV=40');
% xlabel('NS')    

