% Arbtirary example of an American option using the Beliaeva-Nawalkha tree

clc; 
clear;

% Option features
T  = 1/12;
S0 = 100;
rf = 0.05;
q = 0;
PutCall = 'P';
Strike = 100;

% Heston parameters
kappa = 3;
theta = 0.04;
sigma = 0.9;
rho = -0.7;
V0 = 0.06;
trap = 1;
param = [kappa theta sigma V0 rho];

% Tree settings
threshold = 1e-20;
NT = 30;
dt = T/NT;

% Generate the closed European price
[x w] = GenerateGaussLaguerre(32);
EuroClosedPrice = HestonPriceGaussLaguerre(PutCall,S0,Strike,T,rf,q,param,trap,x,w);


%% Generate the bivariate (S,V) tree
% BuildBNTree  uses cells for the branches and probabilities
% BuildBNTree2 uses matrices for the branches, cells for the probabilities
% BuildBNTree3 uses matrices for the branches and probabilities

tic
%[EuroPrice AmerPrice Euro Amer Yt V X Prob Branch] = BuildBivariateTree(S0,PutCall,Strike,T,rf,NT,kappa,theta,sigma,V0,rho,threshold);
%[EuroPrice AmerPrice Euro Amer Yt V X Prob Branch] = BuildBivariateTree2(S0,PutCall,Strike,T,rf,NT,kappa,theta,sigma,V0,rho,threshold);
[EuroPrice AmerPrice Euro Amer Yt V X Prob Branch] = BuildBivariateTree3(S0,PutCall,Strike,T,rf,NT,kappa,theta,sigma,V0,rho,threshold);
toc

% American price by control variate method
CVAmerPrice = EuroClosedPrice + (AmerPrice - EuroPrice);

%% Output the results
fprintf('Heston American option Belieava-Nawalkha tree\n');
fprintf('Using %3.0f time steps \n',NT);
fprintf('----------------------------------\n');
fprintf('Closed European Price %10.6f \n', EuroClosedPrice)
fprintf('Tree   European Price %10.6f \n', EuroPrice)
fprintf('Tree   American Price %10.6f \n', AmerPrice)
fprintf('C.V.   American Price %10.6f \n', CVAmerPrice)
fprintf('----------------------------------\n');



