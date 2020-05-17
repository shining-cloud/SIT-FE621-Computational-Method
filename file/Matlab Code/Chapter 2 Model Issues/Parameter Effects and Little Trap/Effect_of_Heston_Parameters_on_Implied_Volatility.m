% Heston option prices generate the implied volatility smile.
% Using the original formulation by Heston (1993).

clc; clear;

S = 100;         % Spot price.
r = 0.05;        % Risk free rate.
q = 0;           % Dividend yield
tau = .25;       % Time to maturity.
rho =  0.0;      % Heston parameter: correlation.
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.01;    % Heston parameter: mean reversion level.
lambda = 0;      % Heston parameter: risk.
v0 = 0.05;       % Heston parameter: initial variance.

PutCall = 'C';   % 'C'all or 'P'ut.
trap = 1;        % 1 = "Little Trap" form, 0 = Original Heston form
K = [95:0.5:110]; % Strike range

% Settings for the bisection algorithm
a = 0.01;
b = 2;
Tol = 1e-4;
MaxIter = 1000;

% Generate the Gauss Laguerre abscissas and weights
[x w] = GenerateGaussLaguerre(32);

%% Effect of sigma on the IV smile
sigma = [0.35 0.55 0.8];
for j=1:length(sigma);
	for i=1:length(K);
		HCallS(i,j) = HestonPriceGaussLaguerre(PutCall,S,K(i),tau,r,q,kappa,theta,sigma(j),lambda,v0,rho,trap,x,w);
		IVS(i,j) = BisecBSIV(PutCall,S,K(i),r,q,tau,a,b,HCallS(i,j),Tol,MaxIter);
	end
end
	
%% Effect of rho on the IV smile
rho = [-0.25 0 0.25];
sigmaV = 0.25;
for j=1:length(rho);
	for i=1:length(K);
		HCallR(i,j) = HestonPriceGaussLaguerre(PutCall,S,K(i),tau,r,q,kappa,theta,sigmaV,lambda,v0,rho(j),trap,x,w);
		IVR(i,j) = BisecBSIV(PutCall,S,K(i),r,q,tau,a,b,HCallR(i,j),Tol,MaxIter);
	end
end

%% Effect of v0 on the IV smile
v0 = [0.02 0.015 0.01];
sigmaV = 0.5;
rhoV = 0;
for j=1:length(v0);
	for i=1:length(K);
		HCallV(i,j) = HestonPriceGaussLaguerre(PutCall,S,K(i),tau,r,q,kappa,theta,sigmaV,lambda,v0(j),rhoV,trap,x,w);
		IVV(i,j) = BisecBSIV(PutCall,S,K(i),r,q,tau,a,b,HCallV(i,j),Tol,MaxIter);
	end
end

%% Effect of kappa on the IV smile
kappa = [5 2 1];
sigmaV = 0.2;
rhoV = 0.0;
v0V = 0.01;
for j=1:length(kappa);
	for i=1:length(K);
		HCallK(i,j) = HestonPriceGaussLaguerre(PutCall,S,K(i),tau,r,q,kappa(j),theta,sigmaV,lambda,v0V,rhoV,trap,x,w);
		IVK(i,j) = BisecBSIV(PutCall,S,K(i),r,q,tau,a,b,HCallK(i,j),Tol,MaxIter);
	end
end

%% Effect of theta on the IV smile
theta = [0.02 0.015 0.01];
sigmaV = 0.5;
rhoV = 0;
v0V = 0.05;
kappaV = 2;
for j=1:length(theta);
	for i=1:length(K);
		HCallT(i,j) = HestonPriceGaussLaguerre(PutCall,S,K(i),tau,r,q,kappaV,theta(j),sigmaV,lambda,v0V,rhoV,trap,x,w);
		IVT(i,j) = BisecBSIV(PutCall,S,K(i),r,q,tau,a,b,HCallT(i,j),Tol,MaxIter);
	end
end

%% Plot the effect of the parameters

subplot(3,2,1)
plot(K,IVR(:,1),'kx-',K,IVR(:,2),'rx-',K,IVR(:,3),'gx-')
legend(['Rho ' num2str(rho(1),'%5.1f')], ...
       ['Rho ' num2str(rho(2),'%5.1f')], ...
       ['Rho ' num2str(rho(3),'%5.1f')])

% Plot the effect of sigma
subplot(3,2,2)
plot(K,IVS(:,1),'kx-',K,IVS(:,2),'rx-',K,IVS(:,3),'gx-')
legend(['Sigma ' num2str(sigma(1),'%5.2f')], ...
       ['Sigma ' num2str(sigma(2),'%5.2f')], ... 
       ['Sigma ' num2str(sigma(3),'%5.2f')])
ylim([0.18 0.21])
   
% Plot the effect of kappa
subplot(3,2,3)
plot(K,IVK(:,1),'kx-',K,IVK(:,2),'rx-',K,IVK(:,3),'gx-')
legend(['kappa ' num2str(kappa(1),'%5.1f') ],...
       ['kappa ' num2str(kappa(2),'%5.1f') ],...
       ['kappa ' num2str(kappa(3),'%5.1f') ])
   
% Plot the effect of theta
subplot(3,2,4)
plot(K,IVT(:,1),'kx-',K,IVT(:,2),'rx-',K,IVT(:,3),'gx-')
legend(['theta ' num2str(theta(1),'%5.3f')],...
       ['theta ' num2str(theta(2),'%5.3f')],...
       ['theta ' num2str(theta(3),'%5.3f')])
   
% Plot the effect of v0
subplot(3,2,5)
plot(K,IVV(:,1),'kx-',K,IVV(:,2),'rx-',K,IVV(:,3),'gx-')
legend(['v(0) ' num2str(v0(1),'%5.3f')], ...
       ['v(0) ' num2str(v0(2),'%5.3f')], ...
       ['v(0) ' num2str(v0(3),'%5.3f')])

   