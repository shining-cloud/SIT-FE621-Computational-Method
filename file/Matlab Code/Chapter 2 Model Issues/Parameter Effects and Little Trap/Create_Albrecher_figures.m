% Illustration of discontinuities in the Heston integrand
% Reproduces the figures from Albrecher et al. (2006) "The Little Heston Trap"

clc; clear;

% Option features
S = 100;              % Spot price
K = 100;              % Strike price
tau = 5;              % Maturity
r = 0.035;            % Risk free rate
kappa = 1.5768;       % Heston parameter : reversion speed
sigma = 0.5751;       % Heston parameter : volatility of variance
rho   = -0.5711;      % Heston parameter : correlation
theta = 0.0398;       % Heston parameter : reversion level
v0    = 0.0175;       % Heston parameter : initial variance
lambda = 0;           % Heston parameter : risk preference

%% Illustration of the integrand for P1
Pnum = 1;
phi = 0.0001:.1:10;
N = length(phi);
for x=1:N
	HestonP1(x)    = HestonIntegrand(phi(x),kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Pnum,0);
	AlbrecherP1(x) = HestonIntegrand(phi(x),kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Pnum,1);
end;

% Find the points of discontinuity
NN = [];
for x=2:N
    if sign(HestonP1(x)) ~= sign(HestonP1(x-1))
        NN = [NN x-1];
    end
end
N1 = NN(1);
N2 = NN(2);

phi1 = phi(1:N1);
phi2 = phi(N1+1:N);
Part1 = HestonP1(1:N1);
Part2 = HestonP1(N1+1:N);
plot(phi,AlbrecherP1,'r-',phi1,Part1,'k',phi2,Part2,'k-')
h = line([phi(N1) phi(N1)], [HestonP1(N1+1) HestonP1(N1)]);
set(h,'LineStyle','--','Color','k')
legend( 'Albrecher Formulation for P1','Heston Formulation for P1')
xlabel('Integration Range')
ylabel('Integrand')

%% Illustration of the integrand for P2
Pnum = 2;
for x=1:N
	HestonP2(x)    = HestonIntegrand(phi(x),kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Pnum,0);
	AlbrecherP2(x) = HestonIntegrand(phi(x),kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,Pnum,1);
end;

% Find the points of discontinuity
threshold = 0.1;
NN = [];
for x=2:N
    if abs(HestonP2(x) - HestonP2(x-1)) > threshold
        NN = [NN x-1];
    end
end
N1 = NN;

phi1 = phi(1:N1);
phi2 = phi(N1+1:N);
Part1 = HestonP2(1:N1);
Part2 = HestonP2(N1+1:N);
plot(phi,AlbrecherP2,'r-',phi1,Part1,'k',phi2,Part2,'k-')
h = line([phi(N1) phi(N1)], [HestonP2(N1+1) HestonP2(N1)]);
set(h,'LineStyle','--','Color','k')
legend( 'Albrecher Formulation for P2','Heston Formulation for P2')
xlabel('Integration Range')
ylabel('Integrand')

