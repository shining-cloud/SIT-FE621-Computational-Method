% Plot of the Heston integrand
clc; clear;

% Option features
S = 15;          % Spot price
K = 10;         % Strike Price
r = 0.0;        % Risk free rate
q = 0.0;        % Dividend yield
Pnum = 2;       % Characteristic function (j=1,2)
trap = 0;       % Heston formulation for characteristic function

% Heston parameters
lambda = 0;      % Risk preference
kappa = 10;      % Mean reversion speed
theta = 0.07;    % Mean reversion level
sigma = 0.3;     % Volatility of variance
v0    = 0.07;    % Initial variance
rho = -0.9;      % Correlation

% Vector of maturities
T = [0.02:.01:.25];

% Integrand grid
dphi = .01;      % Grid increment
Uphi = 50;      % Grid upper limit
Lphi = -50;     % Grid lower limit
phi = [Lphi:dphi:Uphi];


%% Create the integrand for each maturity
for t = 1:length(T);
    for x = 1:length(phi);
        Int(t,x) = HestonProb(phi(x),kappa,theta,lambda,rho,sigma,T(t),K,S,r,q,v0,Pnum,trap);
    end
end

%% Surface plot for the result
surf(Int,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
axis tight
camlight('left')
camlight left

set(gca,'YTickLabel', [0.01:0.05:.25]);
ylabel('Maturity')
numYticks = 5;
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),numYticks))

set(gca,'XTickLabel', [-50:25:50]);
xlabel('Integrand')
numXticks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),numXticks))

    