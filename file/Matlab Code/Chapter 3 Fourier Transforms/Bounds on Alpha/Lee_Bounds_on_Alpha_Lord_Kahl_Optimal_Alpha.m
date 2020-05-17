% Analysis of Alpha in Carr-Madan Fourier representation
%  Roger Lee Bounds on Alpha
% "Option Pricing by Transform Methods: Extensions, Unification..."
%  Kahl and Lord Optimal alpha 
% "Optimal Fourier Inversion in Semi-Analytical Option Pricing"

clc; clear;

%% Settings for the Heston model
% Uses the settings in Kahl and Lord's paper
S = 1;                   % Spot price
r = 0.0;                 % Risk free rate
K = 1.2;                 % Strike price
tau = 1;                 % Maturity
rho = -0.7;              % Heston parameter: correlation
v0 = 0.1;                % Heston parameter: initial variance
theta = 0.1;             % Heston parameter: mean reversion level
kappa = 1;               % Heston parameter: mean reversion speed
sigma = 1;               % Heston parameter: volatility of variance
lambda = 0;              % Heston parameter: risk preference
trap = 1;                % Characteristic function formulation
                         % 0 = Original Heston formulation
						 % 1 = "Little Trap" formulation (Albrecher et al.)
% Verify : kappa - rho*sigma >0
check = kappa - rho*sigma;

%% Bounds on alpha and optimal alpha 
% Shows multiple solutions to Rogers Lee's formula 
% in which we solve for "a" in g(-ia)*exp(d(-ia)*T) = 1
clear A y z
A = [-15:.01:15];
for j=1:length(A);
	g = RogerLeeG(-i*A(j),kappa,rho,sigma);
	d = RogerLeeD(-i*A(j),kappa,rho,sigma);
	E = real(g*exp(d*tau));
	y(j) = (E-1);
	z(j) = 0;
end

% Plot the multiple solutions
plot(A,y,A,z)
axis([A(1) A(end) -5 35])
xlabel('a')
fprintf('Hit a key to advance to the next graph\n')
pause

%% Find yMax and yMin from Roger Lee's closed form
ymax = (sigma - 2*kappa*rho + sqrt(sigma^2 - 4*kappa*rho*sigma + 4*kappa^2))...
	 / (2*sigma*(1-rho^2));
ymin = (sigma - 2*kappa*rho - sqrt(sigma^2 - 4*kappa*rho*sigma + 4*kappa^2))...
     / (2*sigma*(1-rho^2));

%% Find yMax through optimization of Roger Lee's function
% Look for solutions < ymin and > ymax
% Find y-
start = ymin - 1;
[yneg feval] = fminsearch(@(a) RogerLeeGExpD(-i*a,kappa,rho,sigma,tau), start);

% Find y+
start = ymax + 5;
[ypos feval] = fminsearch(@(a) RogerLeeGExpD(-i*a,kappa,rho,sigma,tau), start);

% The range for Ax
Ax = [yneg ypos];
 
%% Roger Lee bounds: max alpha, min alpha, allowable range for alpha
AlphaMax = ypos - 1;
AlphaMin = yneg - 1;
AlphaRange = Ax - 1;

% Lord and Kahl optimal alpha
start = (AlphaMin + AlphaMax)/2;
[AlphaOptimal eval] = fminsearch(@(alpha) LordKahlFindAlpha(alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap),start);

%% Display the resunts
clc;
fprintf('Roger Lee bounds on alpha\n')
fprintf('-----------------------------------------------------------\n')
fprintf('Ymin and Ymax are          (%7.4f,%7.4f) \n',ymin,ymax )
fprintf('The range Ax is (y-,y+)    (%7.4f,%7.4f) \n',yneg,ypos)
fprintf('(AlphaMin,AlphaMax) is     (%7.4f,%7.4f) = Ax - 1 \n',AlphaMin,AlphaMax)
fprintf('-----------------------------------------------------------\n')
fprintf('Lord and Kahl Optimal alpha %7.4f \n', AlphaOptimal)
fprintf('--------------------------------------------------\n')
fprintf('Note that alpha optimal in (%7.4f,%7.4f) \n',AlphaMin,AlphaMax)
fprintf('-----------------------------------------------------------\n')

%% Plot Lord and Kahl's optimal alpha function and its derivative
% Figure 3(A) of Lord and Kahl
dA = 0.01;
A = [AlphaMin:dA:AlphaMax];
N = length(A);

% The function for optimal alpha
for x=1:N
	f(x) = LordKahlFindAlpha(A(x),kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
end

% The derivative of the function
for x=2:N-1
	der(x) = (f(x+1) - f(x-1))/(2*dA);
end
der = der(2:N-1);

% Find the points where the derivative switches sign.  
% These are the local minima
optim = zeros(N-3,1);
for x=2:N-2
	if sign(der(x)) ~= sign(der(x-1)) & abs(der(x)) < .1;
		optim(x) = 1;
	end
end

% Identify the local minima along the x-axis
Op = find(optim==1);
Opy = zeros(1,length(Op));
Opx = (A(Op) + A(Op+1))/2;
AlphaOptimalSet = Opx;
fstring = repmat('%10.4f',1,length(AlphaOptimalSet));
fprintf(['The optimal alpha set is ' fstring '\n'], AlphaOptimalSet)

% Plot the function, its derivative, and the local minima
z = zeros(length(A),1);

plot(A, f,'k-',A(1:N-2),der,'r-',Opx,Opy,'bo',A,z,'k')
axis([A(1) A(end) -8 5])
legend('The optimal alpha function', 'The derivative of the function', 'Optimal alpha','Location','SouthEast')
fprintf('Hit a key to advance to the next graph\n')
pause

%% Plot Carr-Madan integrand of the Heston model using 3 values for alpha
alpha1 = (AlphaOptimal+AlphaMin)/2;
alpha2 =  AlphaOptimal;
alpha3 = (AlphaOptimal+AlphaMax)/2;

% Figure 2 of Lord and Kahl showing integrand is best behaved
% at the optimal alpha
dv = .01;
V = [.0001:dv:15];
for j=1:length(V);
	v = V(j) - (alpha1+1)*i;
    psi = HestonPsi(V(j),alpha1,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
    psi_1(j) = real(psi);

	v = V(j) - (alpha2+1)*i;
    psi = HestonPsi(V(j),alpha2,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
	psi_2(j) = real(psi);
	
	v = V(j) - (alpha3+1)*i;
    psi = HestonPsi(V(j),alpha3,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
	psi_3(j) = real(psi);
end

plot(V, psi_1,'k-',V,psi_2,'r-',V,psi_3, 'k:');
legend(['Alpha < optimal = ' num2str(alpha1)] , ...
       ['Optimal alpha    = ' num2str(alpha2)] , ...
       ['Alpha > optimal = ' num2str(alpha3)])

