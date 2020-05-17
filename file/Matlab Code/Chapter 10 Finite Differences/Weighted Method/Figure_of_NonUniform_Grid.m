% Plot of a non uniform grid
clc; clear;

% Strike price
K = 10;

% Minimum and maximum values for the Stock Price, Volatility, and Maturity
Smin = 0; Smax = 2.5*K;
Vmin = 0; Vmax = 0.5;

% Number of stock, vol, and time points
nS = 69;
nV = 49;

%% Choose the grid type and build the grid
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
    % The stock price grid
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

%% Plot the grid
plot([0,max(V)],[S',S'],'b')
line([V',V'],[0,max(S)],'Color','r')

