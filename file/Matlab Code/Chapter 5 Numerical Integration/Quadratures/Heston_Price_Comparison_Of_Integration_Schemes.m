% Comparison of Heston (1993) call price obtained by Newton-Coates Formulas
% and obtained by Gauss Legendre and Gauss Laguerre schemes.
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.

clc; clear;

% Option settings
S = 100;         % Spot price.
T = 1/24;        % Time to maturity.
r = 0.03;        % Risk free rate.
q = 0.07;        % Dividend yield
trap = 1;

% Heston parameters
kappa =  8.9832;
theta =  0.0524;
sigma =  1.0982;
v0    =  0.0325;
rho   = -0.9921;
lambda = 0;

% Integration limits
A = 0;           % Lower limit for Gauss Legendre
B = 100;          % Upper limit for Gauss Legendre
a = 1e-10;       % Lower limit for Newton-Cotes
b = 100;          % Upper limit for Newton-Cotes

% Number of points
N = 100;         % Newton-Cotes
NGLa = 5 ;       % Gauss Laguerre
NGLe = 5 ;       % Gauss Legendre
NGLo = 5 ;       % Gauss Lobatto

Method = {'Mid-Point     ' 'Trapezoidal   ' 'Simpsons      ' ...
	      'Simpsons 3/8  ' 'Gauss-Laguerre' 'Gauss-Legendre' ...
		  'Gauss Lobatto '};

% Range of strikes
K = [98:1:102];

% Generate N points for the Gauss Laguerre rule
[xGLa wGLa] =  GenerateGaussLaguerre(NGLa);
[xGLe wGLe] =  GenerateGaussLegendre(NGLe);
[xGLo wGLo] =  GenerateGaussLobatto(NGLe);

%% Calculate the prices and display the results
for k=1:length(K);
    % True prices
    Call(k) = HestonPriceNewtonCoates('C',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,3,1e-5,150,2500);
    Put(k)  = HestonPriceNewtonCoates('P',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,3,1e-5,150,2500);
    % Newton-Cotes formulas
	for method=1:4;
		CallNC(method,k) = HestonPriceNewtonCoates('C',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);
		PutNC(method,k)  = HestonPriceNewtonCoates('P',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N);
        CallError(method,k) = abs(CallNC(method,k) - Call(k));
        PutError(method,k)  = abs(PutNC(method,k) - Put(k));
	end
	% Gauss-Laguerre rules
	CallGLa(k) = HestonPriceGaussLaguerre('C',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLa,wGLa);
	PutGLa(k)  = HestonPriceGaussLaguerre('P',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLa,wGLa);
    CallError(5,k) = abs(CallGLa(k) - Call(k));
    PutError(5,k)  = abs(PutGLa(k) - Put(k));
    % Gauss-Legendre rules
	CallGLe(k) = HestonPriceGaussLegendre('C',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLe,wGLe,A,B);
	PutGLe(k)  = HestonPriceGaussLegendre('P',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLe,wGLe,A,B);
    CallError(6,k) = abs(CallGLe(k) - Call(k));
    PutError(6,k)  = abs(PutGLe(k) - Put(k));
    % Gauss-Lobatto rules
	CallGLo(k) = HestonPriceGaussLegendre('C',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLo,wGLo,a,b);
	PutGLo(k)  = HestonPriceGaussLegendre('P',S,K(k),T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,xGLo,wGLo,a,b);
    CallError(7,k) = abs(CallGLo(k) - Call(k));
    PutError(7,k)  = abs(PutGLo(k) - Put(k));
end

MeanCallError = mean(CallError');
MeanPutError  = mean(PutError');


%% Output the results
fprintf('                   ---------------  Strikes ----------------\n');
fprintf('Call Prices  %10.0f %8.0f %8.0f %8.0f %8.0f    Mean Error\n',K(1),K(2),K(3),K(4),K(5))
fprintf('------------------------------------------------------------------------\n');
fprintf('True Price      %8.2f %8.2f %8.2f %8.2f %8.2f\n',Call(1),Call(2),Call(3),Call(4),Call(5))
fprintf('Mid-Point       %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',CallNC(1,1),CallNC(1,2),CallNC(1,3),CallNC(1,4),CallNC(1,5),MeanCallError(1))
fprintf('Trapezoidal     %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',CallNC(2,1),CallNC(2,2),CallNC(2,3),CallNC(2,4),CallNC(2,5),MeanCallError(2))
fprintf('Simpsons        %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',CallNC(3,1),CallNC(3,2),CallNC(3,3),CallNC(3,4),CallNC(3,5),MeanCallError(3))
fprintf('Simpsons 3/8    %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',CallNC(4,1),CallNC(4,2),CallNC(4,3),CallNC(4,4),CallNC(4,5),MeanCallError(4))
fprintf('Gauss Laguerre  %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',CallGLa(1),CallGLa(2),CallGLa(3),CallGLa(4),CallGLa(5),MeanCallError(5))
fprintf('Gauss Legendre  %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',CallGLe(1),CallGLe(2),CallGLe(3),CallGLe(4),CallGLe(5),MeanCallError(6))
fprintf('Gauss Lobatto   %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',CallGLo(1),CallGLo(2),CallGLo(3),CallGLo(4),CallGLo(5),MeanCallError(7))
fprintf('------------------------------------------------------------------------\n');
fprintf(' \n');

fprintf('                   ---------------  Strikes ----------------\n');
fprintf('Put Prices   %10.0f %8.0f %8.0f %8.0f %8.0f    Mean Error\n',K(1),K(2),K(3),K(4),K(5))
fprintf('------------------------------------------------------------------------\n');
fprintf('True Price      %8.2f %8.2f %8.2f %8.2f %8.2f \n',Put(1),Put(2),Put(3),Put(4),Put(5))
fprintf('Mid-Point       %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',PutNC(1,1),PutNC(1,2),PutNC(1,3),PutNC(1,4),PutNC(1,5),MeanPutError(1))
fprintf('Trapezoidal     %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',PutNC(2,1),PutNC(2,2),PutNC(2,3),PutNC(2,4),PutNC(2,5),MeanPutError(2))
fprintf('Simpsons        %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',PutNC(3,1),PutNC(3,2),PutNC(3,3),PutNC(3,4),PutNC(3,5),MeanPutError(3))
fprintf('Simpsons 3/8    %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',PutNC(4,1),PutNC(4,2),PutNC(4,3),PutNC(4,4),PutNC(4,5),MeanPutError(4))
fprintf('Gauss Laguerre  %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',PutGLa(1),PutGLa(2),PutGLa(3),PutGLa(4),PutGLa(5),MeanPutError(5))
fprintf('Gauss Legendre  %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',PutGLe(1),PutGLe(2),PutGLe(3),PutGLe(4),PutGLe(5),MeanPutError(6))
fprintf('Gauss Lobatto   %8.2f %8.2f %8.2f %8.2f %8.2f %10.4f\n',PutGLo(1),PutGLo(2),PutGLo(3),PutGLo(4),PutGLo(5),MeanPutError(7))
fprintf('------------------------------------------------------------------------\n');


