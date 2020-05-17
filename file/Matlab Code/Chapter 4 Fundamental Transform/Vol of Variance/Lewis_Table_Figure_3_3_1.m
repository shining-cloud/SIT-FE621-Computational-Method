% Volatility of volatility expansion for the Heston call price
% from Alan Lewis "Option Valuation Under Stochastic Volatility" (2001)
% Table 3.3.1 of Lewis (2001) page 81 and Figure 3.3.1. on page 80.

clc; clear;

%% Parameter and other settings
% Input the Heston parameters.  These are identical to those in table 3.3.1
% on page 81 of Lewis (2001) but reparameterized as in Heston (1993)
lambda = 0;         % Risk  parameter
rho = -0.5;         % Correlation
v0 = 0.0225;        % Initial variance
kappa = 4;          % Mean reversion speed
theta = 0.09/4;     % Mean reversion level
sigma = 0.1;        % Volatility of variance
 
% Option features
T  = 0.25;           % Time to maturity in years
rf = 0.0;            % Risk free rate
q  = 0.0;            % Dividend yield
S  = 100;            % Spot price
K  = 100;            % Strike price
trap = 1;            % Formulation of the characteristic function
                     % 1 = "Little Heston Trap" form
                     % 2 =  original Heston formulation

% Specify Newton Coates method
%   1 = mid-point rule, 2 = trapezoidal rule
%   3 = Simpson's rule, 4 = Simpson's 3/8 rule
method = 4;

%% Find the Heston prices using exact method and series expansions
% Find the exact Heston call price
a = 1e-10;
b = 100;
N = 10000;

BSC = @(S,K,rf,q,v,T) (S*exp(-q*T)*normcdf((log(S/K) + T*(rf - q + 0.5*v)) / sqrt(v*T)) - exp(-rf*T)*K*normcdf((log(S/K) + T*(rf - q + 0.5*v)) / sqrt(v*T) - sqrt(v*T)));
		 
% Time average of the deterministic variance
v = theta + (v0-theta)*(1-exp(-kappa*T))/(kappa*T);

% Time average of the deterministic volatility
IV = sqrt(v);

%% Reproduce Table 3.3.1 on Page 81 of Lewis' book (2001)
K = [70:10:130];
T  = 0.25;
A = 0.001;
B = 10;
Tol = 1e-5;
MaxIter = 1000;

for k=1:length(K)
    % Find the implied vol for the Series I price using the BSIV function
    SeriesICallPrice(k) = SeriesICall(S,K(k),rf,q,T,v0,rho,theta,kappa,sigma);
    IV1(k) = BisecBSIV('C',S,K(k),rf,q,T,A,B,SeriesICallPrice(k),Tol,MaxIter);
    % Find the implied volatility for the Series II price directly
    [IV2(k) SeriesIICallPrice(k)] = SeriesIICall(S,K(k),rf,q,T,v0,rho,theta,kappa,sigma);
    Exact(k) = HestonPriceNewtonCoates('C',S,K(k),T,rf,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);
    IVExact(k) = BisecBSIV('C',S,K(k),rf,q,T,A,B,Exact(k),Tol,MaxIter);
    % Black Scholes price
    BSCall(k) = BSC(S,K(k),rf,q,v,T);
end

%% Display the results
clc;
disp(['Strike Price            ' num2str(K, '%10.0f')])
disp('---------------------------------------------------------------------------------------------')
disp(['Exact Price             ' num2str(Exact,'%10.6f')])
disp(['Implied Volatility      ' num2str(100.*IVExact,'%10.2f')])
disp(['Series I CallPrice      ' num2str(SeriesICallPrice, '%10.6f')])
disp(['Series I Implied Vol    ' num2str(100.*IV1, '%10.2f')])
disp(['Series II CallPrice     ' num2str(SeriesIICallPrice, '%10.6f')])
disp(['Series II Implied Vol   ' num2str(100.*IV2, '%10.2f')])
disp(['Black Scholes Price     ' num2str(BSCall, '%10.6f')])
disp(['Black Scholes Imp Vol   ' num2str(100.*IV.*ones(1,7),'%10.2f')])
disp('---------------------------------------------------------------------------------------------')


%% Reproduce Figure 3.3.1 on Page 80 of Lewis' book (2001)
clear K IV1 IV2 SeriesICallPrice Exact IVExact
T = [0.25 1.5];
K = [75:125]';

for t=1:length(T)
	for k=1:length(K)
		% Find the implied vol for the Series I price using the BSIV function
		SeriesICallPrice(k,t) = SeriesICall(S,K(k),rf,q,T(t),v0,rho,theta,kappa,sigma);
		IV1(k,t) = BisecBSIV('C',S,K(k),rf,q,T(t),A,B,SeriesICallPrice(k,t),Tol,MaxIter);
		% Find the implied volatility for the Series II price directly
		[IV2(k,t) junk] = SeriesIICall(S,K(k),rf,q,T(t),v0,rho,theta,kappa,sigma);
		% Find the implied volatility from the exact price
		Exact(k,t) = HestonPriceNewtonCoates('C',S,K(k),T(t),rf,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);
		IVExact(k,t) = BisecBSIV('C',S,K(k),rf,q,T(t),A,B,Exact(k,t),Tol,MaxIter);
	end
end

%% Plot the results
plot(K,IV1(:,1),'r-',...
     K,IV2(:,1),'b-',...
     K,IVExact(:,1),'k-',...
     K,IV1(:,2),'rx-',...
     K,IV2(:,2),'bx-',...
     K,IVExact(:,2),'kx-')
legend('Series I  Implied Vol,  3 months mat.',...
       'Series II Implied Vol,  3 months mat.',...
       'Exact Price Implied Vol, 3 months mat.',...
       'Series I  Implied Vol, 18 months mat.',...
       'Series II Implied Vol, 18 months mat.',...
       'Exact Price Implied Vol, 18 months mat.');
xlabel('Strike Price')
ylabel('Implied Volatility')
