% Plot of the Heston integrand

clc; clear;

K = 75;         % Strike Price
S = 100;        % Spot Price
r = 0.0;        % Risk free rate
q = 0.0;        % Dividend yield
Pnum = 1;       % Characteristic function (j=1,2)
Trap = 0;       % 1=Little Trap form, 0=Heston form

kappa = 10;     % Heston parameter: mean reversion speed
theta = 0.05;   % Heston parameter: mean reversion level
v0 = 0.05;      % Heston parameter: initial variance
rho = -0.9;     % Heston parameter: correlation
lambda = 0;     % Heston parameter: risk preference

sigma1 = 0.75;  % Heston parameter: volatility of variance
tau1 = 3;       % Time to maturity

sigma2 = 0.09;  % Heston parameter: volatility of variance
tau2 = 1;       % Time to maturity
Trap2 = 1;

dphi = 0.001;   % Integration grid increment
Uphi = 10;      % Integration grid upper limit
Lphi = 0.00001; % Integration grid lower limit

PHI = [Lphi:dphi:Uphi];
N = length(PHI);


%% Obtain the integrals with the Heston form of the characteristic function
for x=1:N
	phi = PHI(x);
	Inte1(x) = HestonProb(phi,kappa,theta,lambda,rho,sigma1,tau1,K,S,r,q,v0,Pnum,Trap);
	Inte2(x) = HestonProb(phi,kappa,theta,lambda,rho,sigma2,tau2,K,S,r,q,v0,Pnum,Trap);
end


%% Plot the results, showing the discontinuity

% Points of discontinuity
N1 = 1623;
N2 = 1624;
N3 = 5094;
N4 = 5095;

% Different portions of the integration range and integrand
phi1 = PHI(1:N1);
phi2 = PHI(N2:N3);
phi3 = PHI(N4:N);
Part1 = Inte1(1:N1);
Part2 = Inte1(N2:N3);
Part3 = Inte1(N4:N);

% The plot
plot(PHI,Inte2,'r-',phi1,Part1,'k-',phi2,Part2,'k-',phi3,Part3,'k-');
h1 = line([PHI(N1) PHI(N1)], [Inte1(N1) Inte1(N2)]);
h2 = line([PHI(N3) PHI(N3)], [Inte1(N3) Inte1(N4)]);
set(h1,'LineStyle','--','Color','k')
set(h2,'LineStyle','--','Color','k')
legend(['Maturity ' num2str(tau1) ' years, sigma = ' num2str(sigma1)], ...
	   ['Maturity ' num2str(tau2) ' year, sigma = ' num2str(sigma2)])
xlabel('Integration range')
ylabel('Integrand')

