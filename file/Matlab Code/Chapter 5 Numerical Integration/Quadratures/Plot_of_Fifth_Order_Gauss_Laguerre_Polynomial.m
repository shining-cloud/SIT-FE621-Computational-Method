% Plot the Laguerre Polynomial of order 5
clc; clear;

% Laguerre polynomial coefficients and roots
C = [-1/120 5/24 -5/3 5 -5 1];
R  = flipud(roots(C));
degree = length(C)-1;

% Domain
X = [0:0.01:12.8];
N = length(X);
X0 = [0:14];
N0 = length(X0);

%% Evaluate the polynomial at each point in the domain
for x=1:N
    L(x) = polyval(C,X(x));
end

% Plot the result
Z0 = zeros(1,N0);
ZR = zeros(1,length(R));
plot(X,L,'r-',R,ZR,'ko',X0,Z0,'k-')
legend(['Laguerre polynomial of degree ' num2str(degree)],'Polynomial roots');
   
    