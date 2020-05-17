% Alternating Direction Schemes
% Examples of schemes
clc; clear;

% Set the scheme 'DO' Douglas, 'CS' Craig-Sneyd,
% 'MCS' Modified Craig-Sneyd, 'HV' Hundsdorfer-Verview
scheme = {'DO' 'CS' 'MCS' 'HV'};


% Strike price, risk free rate, dividend yield, and maturity
K   = 100;
r   = 0.02;
q   = 0.05;
Mat = 0.15;

% Heston parameters.  Case 1 of Hout and Foulon (Table 1)
S0 = 101.52;
V0 = 0.05412;
kappa =  1.5;
theta =  0.04;
sigma =  0.3;
rho   = -0.9;
v0    =  0.05;
lambda = 0;
params = [kappa theta sigma V0 rho lambda];

% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0; Smax = 2*K;
Vmin = 0; Vmax = 0.5;
Tmin = 0; Tmax = Mat;

% Number of stock, volatility, and time points
nS = 19;
nV = 19;
nT = 9;

% Choose the grid type
GridType = 'NonUniform';

if strcmp(GridType,'Uniform')
    % Increment for Stock Price, Volatility, and Maturity
    ds = (Smax-Smin)/nS;
    dv = (Vmax-Vmin)/nV;
    % Grid Vectors for the Stock Price, Volatility, and Maturity
    NS = nS+1;
    NV = nV+1;
    S = [0:NS-1].*ds;
    V = [0:NV-1].*dv;
else
    % The non uniform stock price grid
    c = K/5;
    dz = 1/nS*(asinh((Smax-K)/c) - asinh(-K/c));
    for i=1:nS+1;
        z(i) = asinh(-K/c) + (i-1)*dz;
        S(i) = K + c*sinh(z(i));
    end
    % The volatility grid
    d = Vmax/10;
    dn = asinh(Vmax/d)/nV;
    for j=1:nV+1
        n(j) = (j-1)*dn;
        V(j) = d*sinh(n(j));
    end
    NS = nS+1;
    NV = nV+1;
end

% Time grid (common to both uniform and non uniform grid types)
NT = nT+1;
dt = (Tmax-Tmin)/nT;
T = [0:NT-1].*dt;


%% Run the ADI schemes

for s=1:length(scheme)
    % Explicit Scheme
    thet = 0;
    disp('Running Explicit ADI Scheme --------')
    EPrice(s) = ADIPrice(scheme(s),thet,params,S0,V0,K,r,q,S,V,T,GridType);
    disp('Done')

    % Implicit Scheme
    thet = 1;
    disp('Running Implicit ADI Scheme --------')
    IPrice(s) = ADIPrice(scheme(s),thet,params,S0,V0,K,r,q,S,V,T,GridType);
    disp('Done')

    % Crank-Nicolson Scheme
    thet = 0.5;
    disp('Running Crank-Nicolson ADI Scheme --------')
    CPrice(s) = ADIPrice(scheme(s),thet,params,S0,V0,K,r,q,S,V,T,GridType);
    disp('Done')
end

%% Calculate the exact price
trap = 1;
[x w] = GenerateGaussLaguerre(32);
HPrice = HestonPriceGaussLaguerre('C',S0,K,Mat,r,q,kappa,theta,sigma,lambda,V0,rho,trap,x,w);

%% Errors
for s=1:length(scheme)
    ExError(s) = (EPrice(s) - HPrice);
    ImError(s) = (IPrice(s) - HPrice);
    CError(s)  = (CPrice(s) - HPrice);
end

%% Output everything
clc;
fprintf('Grid Type %s \n',GridType);
fprintf('Stock price grid size of %2.0f \n', NS);
fprintf('Volatility grid size of  %2.0f \n', NV);
fprintf('Number of time steps     %2.0f \n', NT);
fprintf('Exact Price is           %2.4f \n', HPrice);
fprintf('----------------------------------------------------------------------------------\n')
fprintf('ADI     Explicit    Error     Implicit    Error   Crank-Nicol   Error\n')
fprintf('----------------------------------------------------------------------------------\n')
scheme = {'DO ' 'CS ' 'MCS' 'HV '};
for s=1:length(scheme)
    fprintf('%s  %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',...
            char(scheme(s)),EPrice(s),ExError(s),IPrice(s),ImError(s),CPrice(s),CError(s));
end
fprintf('----------------------------------------------------------------------------------\n')

