% Heston American put prices of Clarke and Parrott (1999)
% Uses Medvedev-Scaillet Expansion

clc; clear;

% Settings from Clarke and Parrott
S = 10;
r = 0.1;
q = 0;

% Heston parameters
kappaV = 5;
thetaV = 0.16;
sigmaV = 0.9;
rho = 0.1;
v0  = 0.0625;
lambda = 0;
trap = 1;
params = [kappaV thetaV sigmaV v0 rho];
PutCall = 'P';


%% Settings for Simpson's Rule
A = 1e-20;
B = 100;
N = 10000;
method = 3;

yinf = 1e4;        % Infinite barrier for the European put
hi = 4;            % Upper bound for optimization

% Select the number of terms in the expansion (3,4 or 5)
NumTerms = 3;

%% Find the finite difference prices and Greeks
T = [0.25 0.5 0.75 1];
K = [8 9 10 11 12];

GreekChoice = 'vega1';

for t=1:length(T);
    for k=1:length(K);
        Greek(k,t) = MSGreeksFD(params,S,K(k),r,q,T(t),method,A,B,N,yinf,hi,NumTerms,GreekChoice);
    end
end

%% Plot the Greeks
% Interpolate along a finer grid
Ki = [K(1):.1:K(end)];
Ti = [T(1):.025:T(end)];
GreekI  = interp2(T,K,Greek, Ti,Ki','spline');

% Axis ticks
Kind = [find(Ki==K(1)) find(Ki==K(2)) find(Ki==K(3)) find(Ki==K(4)) find(Ki==K(5))];
Tind = [find(Ti==T(1)) find(Ti==T(2)) find(Ti==T(3)) find(Ti==T(4)) ];

% Produce the surface plot
mesh(GreekI,'FaceColor','interp','EdgeColor','k','FaceLighting','phong')
xlabel('Maturity')
ylabel('Strike Price')
set(gca,'XTick',      Tind);  % Maturity
set(gca,'XTickLabel', T)        
set(gca,'YTick',      Kind);  % Strike
set(gca,'YTickLabel', K);
axis tight
view([-300,150,70])
set(gca,'YDir','reverse')


