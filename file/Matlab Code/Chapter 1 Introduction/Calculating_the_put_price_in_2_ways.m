% Obtaining the price of the Heston put in two ways:
% By put-call parity, and directly

clc; clear;

% Option features
S = 100;         % Spot price
K = 95;          % Strike price
T = .5;          % Maturity
r = 0.04;        % Risk free rate
q = 0.03;        % Dividend yield
kappa = 5;       % Heston parameter : reversion speed
sigma = 0.5;     % Heston parameter : volatility of variance
rho   = -0.8;    % Heston parameter : correlation
theta = 0.05;    % Heston parameter : reversion level
v0    = 0.05;    % Heston parameter : initial variance
lambda = 0;      % Heston parameter : risk preference
trap = 0;        % 0 = Original Heston formulation
                 % 1 = Albrecher et al formulation

% Integration range				 
Lphi = 1e-10;     % Lower limit
dphi = 0.01;      % Increment
Uphi = 100;       % Upper limit

% Build the integration grid
phi = [Lphi:dphi:Uphi];
N = length(phi);

%% Build the integrands for P1 and P2;
for k=1:N;
	int1(k) = HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
	int2(k) = HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
end

% The integrals
I1 = trapz(int1)*dphi;
I2 = trapz(int2)*dphi;

%% Put price using call price and put call parity
% The probabilities P1 and P2
P1 = 1/2 + 1/pi*I1;
P2 = 1/2 + 1/pi*I2;

% The call price
HestonC = S*exp(-q*T)*P1 - K*exp(-r*T)*P2;

% The put price by put-call parity
HestonPut1 = HestonC - S*exp(-q*T) + K*exp(-r*T);

%% Put price directly
% Redefine the probabilities P1 and P2
p1 = 1/2 - 1/pi*I1;
p2 = 1/2 - 1/pi*I2;

% Put price directly
HestonPut2 = K*exp(-r*T)*p2 - S*exp(-q*T)*p1;

%% Output the results
fprintf('--------------------------------------------------------- \n');
fprintf('Calculation method         P1          P2       PutPrice \n');
fprintf('--------------------------------------------------------- \n');
fprintf('Put-Call Parity       %10.4f  %10.4f  %10.4f \n',P1,P2,HestonPut1);
fprintf('Direct calculation    %10.4f  %10.4f  %10.4f \n',p1,p2,HestonPut2);
fprintf('--------------------------------------------------------- \n');


