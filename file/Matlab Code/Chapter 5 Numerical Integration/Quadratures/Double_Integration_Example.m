% Illustration of the double trapezoidal rule on bivariate normal cdf
clc; clear;

% Parameters of the bivariate normal cdf
mu  = [0 0];
cov = [1 0; 0 1];

% Point at which to evaluate the CDF
Z = [-0.515,0.243];

% True value of the CDF returned by Matlab
TrueValue = mvncdf(Z,mu,cov);

%% Create the bivariate pdf
f = @(x,y) (0.5/pi*exp(-(x^2+y^2)/2));
x = Z(1);
y = Z(2);

% Lower limits for (X,Y)
LoX = -5;
LoY = -5;

% Number of points for the double trapz rule and grid increments
NX = 100;
NY = 100;
dx = (x - LoX)/NX;
dy = (y - LoY)/NY;

% Create the (X,Y) grid
for j=1:NX;
    X(j) = LoX + j*dx;
end
for j=1:NY;
    Y(j) = LoY + j*dy;
end

%% The CDF by the double trapezoidal rule
TrapValue = DoubleTrapz(f,X,Y);


%% The CDF by the double Gauss Legendre rule
[x1 w1] = GenerateGaussLegendre(12);
[x2 w2] = GenerateGaussLegendre(14);
a = LoX;
b = Z(1);
c = LoY;
d = Z(2);
GLeValue = DoubleGaussLegendre(f,a,b,c,d,x1,w1,x2,w2);


%% Output the results
fprintf('Trapezoidal rule using the standard normal bivariate CDF \n');
fprintf('---------------------------------------------------------\n')
fprintf('Matlab value         %10.8f \n',TrueValue);
fprintf('Trapezoidal value    %10.8f \n',TrapValue);
fprintf('Gauss Legendre value %10.8f \n',GLeValue);

