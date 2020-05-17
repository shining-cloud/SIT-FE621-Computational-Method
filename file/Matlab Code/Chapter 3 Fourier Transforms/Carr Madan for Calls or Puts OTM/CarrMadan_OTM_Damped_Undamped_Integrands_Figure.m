% Illustrates integrands for the Carr-Madan representation for OTM options

clc; clear;

% Define the parameters and inputs
S = 1;           % Spot price.
K = 0.96;        % Strike Price
tau = 1/52;      % Time to maturity.
r = 0.03;        % Risk free rate.
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.25;    % Heston parameter: mean reversion level.
sigma = 0.3;     % Heston parameter: volatility of vol
lambda = 0;      % Heston parameter: risk.
v0 = .05;        % Heston parameter: initial variance.
rho = -0.8;      % Heston parameter: correlation
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
alpha = 1.1;     % Carr-Madam dampening factor (for puts, >1)


%% Define two integrands - undamped and damped
U = [-170:170];
T = [1 2 3 4]./52;
for t=1:length(T)
    for u=1:length(U);
        % Undamped integrand
        Undamped(t,u) = CarrMadanIntegrandOTM(U(u),kappa,theta,lambda,rho,sigma,T(t),K,S,r,v0,trap);
        % Damped integrand
        a = CarrMadanIntegrandOTM(U(u)-i*alpha,kappa,theta,lambda,rho,sigma,T(t),K,S,r,v0,trap);
        b = CarrMadanIntegrandOTM(U(u)+i*alpha,kappa,theta,lambda,rho,sigma,T(t),K,S,r,v0,trap);
        Damped(t,u) = (a - b)/2;
    end
end

clear alpha
%% Surface plot

mesh(Damped,'FaceColor','interp','EdgeColor','g','FaceLighting','phong')
hold on
mesh(Undamped,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
alpha(0.5)
hold off
axis tight
ylabel('Maturity (weeks)')        
xlabel('Integration Domain')
set(gca,'XTick',      [1  85 170 255 340]);          % Phi
set(gca,'XTickLabel', [-170 -85 0 85 170])        
set(gca,'YTick',      [1 2 3 4]); % Maturiy
set(gca,'YTickLabel', [4 3 2 1]);
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')


%% Define two integrands - undamped and damped
% U = [-200:.1:200];
% for u=1:length(U);
% 	% Undamped integrand
% 	Undamped(u) = CarrMadanIntegrandOTM(U(u),kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
% 	% Damped integrand 
% 	a(u) = CarrMadanIntegrandOTM(U(u)-i*alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
% 	b(u) = CarrMadanIntegrandOTM(U(u)+i*alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
% 	Damped(u) = (a(u) - b(u))/2;
% end
% 
% % Plot the results
% plot(U,Undamped,'r-',U,Damped,'k-')
% legend('Original integrand', 'Damped Integrand')
