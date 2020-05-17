% Heston American put prices of Clarke and Parrott (1999)
% Uses Medvedev-Scaillet Expansion

clc; clear;

% Settings from Clarke and Parrott
S = [8 9 10 11 12];
TruePrice = [2.00 1.107641 0.520030 0.213668 0.082036];
Strike = 10;
r = 0.1;
q = 0;
T = 0.25;

% Heston parameters
kappaV = 5;
thetaV = 0.16;
sigmaV = 0.9;
rho = 0.1;
v0  = 0.0625;
lambda = 0;
trap = 1;
params = [kappaV thetaV sigmaV v0 rho];

% Settings for Simpson's Rule
A = 1e-20;
B = 100;
N = 1000;
method = 3;

% Infinite barrier for the European put
yinf = 1e4;

% Select the number of terms in the expansion (3,4 or 5)
NumTerms = 4;

%% Find the Medvedev-Scaillet Heston price
for k=1:5
    [EuroPutClosed AmerPutMS AmerPut(k) EEP theta(k) y(k)] = MSPrice(S(k),Strike,T,r,q,params,trap,method,A,B,N,NumTerms,yinf);
end

%% Generate barrier levels for various values of y
Y = [.25:.005:4];
NY = length(Y);
for k=1:5
    for j=1:NY
        B(j,k) = MSPutHeston(Y(j),theta(k),Strike,params,r,q,T,NumTerms);
    end
end
B(:,2) = B(:,2)+ mean(B(:,1) - B(:,2));
B(:,3) = B(:,3)+ mean(B(:,1) - B(:,3));
B(:,4) = B(:,4)+ mean(B(:,1) - B(:,4));
B(:,5) = B(:,5)+ mean(B(:,1) - B(:,5));


%% Plot the barriers 1,2 and 4
plot(Y,B(:,1),'b-',Y,B(:,2),'k-',Y,B(:,4),'r-')
h1 = line([y(1) y(1)],[0 max(B(:,1))]);
h2 = line([y(2) y(2)],[0 max(B(:,2))]);
h4 = line([y(4) y(4)],[0 max(B(:,4))]);
set(h1,'LineStyle',':','Color','b')
set(h2,'LineStyle',':','Color','k')
set(h4,'LineStyle',':','Color','r')
legend(['Spot = ' num2str(S(1)) ',   y = ' num2str(y(1),'%3.2f')],...
       ['Spot = ' num2str(S(2)) ',   y = ' num2str(y(2),'%3.2f')],... 
       ['Spot = ' num2str(S(4)) ', y = '   num2str(y(4),'%3.2f')],...
        'Location','NorthWest')
axis([0 4.25 1.7 2.1])

