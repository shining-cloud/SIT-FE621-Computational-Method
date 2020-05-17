% Variances to reproduce Figure 1 of Chiarella and Ziogas (2006)

clc; clear;

% Expected value of variance
EV = @(Vs,theta,kappa,t,s) (theta + (Vs - theta)*exp(-kappa*(t-s)));

% Maturity and number of time steps
tau = 5;
N = 2500;
dtau = tau/N;

% Heston parameters
kappa = 2;
sigma = 0.1;
theta = 0.01;
v0t = 0.015;

% Initialize variances and maturities
v1(1) = 0.015;
v0(1) = 0.005;
mat(1) = 0;

%% Obtain the values of v0(n) and v1(n)
for n=2:N
    % Maturity
    mat(n) = n*dtau;
    % Variances
    Evt(n) = EV(v0t,theta,kappa,tau,mat(n));
    v0(n) = Evt(n) + sigma/kappa*sqrt(kappa*theta/2);
    v1(n) = Evt(n) - sigma/kappa*sqrt(kappa*theta/2);
end

%% Plot the results, Figure 1 of Chiarella and Ziogas (2006)
X = 2:N;
plot(mat(X),v0(X),'r-',mat(X),Evt(X),'k--',mat(X),v1(X),'b-')
legend('Upper variance v1(t)',...
       'Expected variance E[v(t)]',...
       'Lower variance v0(t)',...
       'Location','NorthWest')
xlabel('Maturity')
ylabel('Variance')
axis([0 5 0 0.025])

