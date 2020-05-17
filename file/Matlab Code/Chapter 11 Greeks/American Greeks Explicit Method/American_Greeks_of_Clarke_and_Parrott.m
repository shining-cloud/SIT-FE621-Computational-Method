% Illustration of Greeks by PDE

clc; clear;

% Import the stock prices and true values
S0 = [8 9 10 11 12];
TruePrice = [2.0000, 1.107641, 0.520030, 0.213668, 0.082036];

% Settings from Clarke and Parrot (1999)
K = 10;
kappa = 5;
theta = 0.16;
sigma = 0.9;
v0 = 0.0625;
rho = 0.1;
lambda = 0;
Mat = 1/4;
r = 0.1;
q = 0.0;
lambda = 0;
params = [kappa theta sigma v0 rho lambda];

% Strike price, risk free rate, dividend yield, and maturity
PutCall = 'P';
EuroAmer = 'A';

% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0;  Smax = 2*K;
Vmin = 0;  Vmax = 0.5;
Tmin = 0;  Tmax = Mat;

% Number of grid points for the stock, volatility, and maturity
nS = 19;        % Stock price
nV = 19;        % Volatility
nT = 3000;      % Maturity

% The maturity time increment and grid
dt = (Tmax-Tmin)/nT;
T = [0:nT].*dt;


%% Obtain the PDE prices and Greeks using a non-uniform grid
for s=1:length(S0);
    % The stock price grid
    c = S0(s)/5;    % Instead of K/5
    dz = 1/nS*(asinh((Smax-S0(s))/c) - asinh(-S0(s)/c));  % Instead of K/c
    for i=1:nS+1;
        z(i) = asinh(-S0(s)/c) + (i-1)*dz;   % Instead of K/c
        S(i) = S0(s) + c*sinh(z(i));
    end

    % The volatility grid
    d = Vmax/500;
    dn = asinh(Vmax/d)/nV;
    for j=1:nV+1
        n(j) = (j-1)*dn;
        V(j) = d*sinh(n(j));
    end

    % Solve the PDE
    [U u] = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer);
    clear dv ds
    dS = 0.01*S0(s);
    dv = 0.01*v0;
    % Price, Delta, and Gamma
    PDEPrice(s) = interp2(V,S,U,v0,S0(s));
    D1 = interp2(V,S,U,v0,S0(s)+dS);
    D2 = interp2(V,S,U,v0,S0(s)-dS);
    DeltaPDE(s) = (D1 - D2)/2/dS;
    GammaPDE(s) = (D1 - 2*PDEPrice(s) + D2)/dS^2;

    % Vega #1
    V1 = interp2(V,S,U,v0+dv,S0(s));
    V2 = interp2(V,S,U,v0-dv,S0(s));
    dCdv0 = (V1 - V2)/2/dv;
    Vega1PDE(s) = dCdv0 * 2.0 * sqrt(v0);
    % Vanna and Volga
    C1 = interp2(V,S,U,v0+dv,S0(s)+dS);
    C2 = interp2(V,S,U,v0-dv,S0(s)+dS);
    C3 = interp2(V,S,U,v0+dv,S0(s)-dS);
    C4 = interp2(V,S,U,v0-dv,S0(s)-dS);
    VannaPDE(s) = (C1 - C2 - C3 + C4)/4.0/dv/dS*2.0*sqrt(v0);
    dC2 = (V1 - 2*PDEPrice(s) + V2)/dv/dv;
    VolgaPDE(s) = 4.0*sqrt(v0)*(dC2*sqrt(v0) + Vega1PDE(s)/4/v0);
    % Theta
    T1 = interp2(V,S,U,v0,S0(s));  % U(S,v,T)
    T2 = interp2(V,S,u,v0,S0(s));  % U(S,v,T-dt)
    Theta(s) = -(T1 - T2)/dt;
end


%% Output the results
clc;
fprintf('Stock price grid size  %5.0f\n', nS+1)
fprintf('Volatility grid size   %5.0f\n', nV+1)
fprintf('Number of time steps   %5.0f\n', nT)
fprintf('---------------------------------------------------------------------------------------\n')
fprintf('Spot TruePrice    PDEPrice    Delta      Gamma      Vega1      Vanna      Volga      Theta\n')
fprintf('---------------------------------------------------------------------------------------\n')
for s=1:length(S0)
    fprintf('%3.0f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n',...
       S0(s),TruePrice(s),PDEPrice(s),DeltaPDE(s),GammaPDE(s),Vega1PDE(s),VannaPDE(s),VolgaPDE(s),Theta(s));
end
fprintf('---------------------------------------------------------------------------------------\n')


%% Validate theta by finite differences
clear Tmax
dT = 0.01*Mat;
Tmax = Mat + dT;

% The maturity time increment and grid
dt = (Tmax-Tmin)/nT;
T = [0:nT].*dt;
[Up u] = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer);

clear Tmax
Tmax = Mat - dT;
dt = (Tmax-Tmin)/nT;
T = [0:nT].*dt;
[Um u] = HestonExplicitPDENonUniformGrid(params,K,r,q,S,V,T,PutCall,EuroAmer);

for s=1:length(S0);
    % Price, Delta, and Gamma
    Pricep(s) = interp2(V,S,Up,v0,S0(s));
    Pricem(s) = interp2(V,S,Um,v0,S0(s));
    ThetaFD(s) = -(Pricep(s) - Pricem(s))/2/dT;
end

fprintf('Theta by F.D.\n')
fprintf('-------------\n')
disp(ThetaFD')
fprintf('-------------\n')

