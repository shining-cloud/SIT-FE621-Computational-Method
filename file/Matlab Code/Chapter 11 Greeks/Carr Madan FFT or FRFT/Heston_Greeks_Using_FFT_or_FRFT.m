% Carr and Madan Greeks using the Fast Fourier Transform
% or the Fractional FFT

clc; clear;

% Required inputs
N = 2^9;         % Number points for the FFT or FRFT
S0 = 100;        % Spot price.
tau = 0.5;       % Time to maturity.
r = 0.05;        % Risk free rate.
q = 0.0;         % Dividend yield
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.06;    % Heston parameter: mean reversion level.
sigma = 0.1;     % Heston parameter: volatility of vol
v0    = 0.06;    % Heston parameter: initial variance.
rho   = -0.7;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
alpha = 1.75;    % Carr-Madam dampening factor
uplimit = 100;   % Upper limit of integration
fast = 1;        % Choice of fast (1) or slow (0) algorithm

%% Run the FFT or FRFT
rule = 'S';
Greek = {'Price','Delta','Gamma','Theta','Rho','Vega1','Vanna','Volga'};
GreeksFFT = zeros(N,8);
eta = 0.1;
lambdainc = .005 ;

% Alternative: Select lambdainc so that the strike range is (Spot +/ 5);
%lambdainc = 2/N*log(S0/(S0-5));

method = 'FRFT';
% FFT or FRFT
for j=1:8;
    if strcmp(method,'FFT')
        [GreeksFFT(:,j) K lambdainc eta] = HestonCallFFTGreek(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,trap,alpha,fast,rule,Greek(j));
    else
        [GreeksFFT(:,j) K lambdainc eta] = HestonCallFRFTGreek(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,trap,alpha,rule,eta,lambdainc,Greek(j));
    end
end

%% Obtain the results near the ATM strikes
ATM = [find(round(K)==S0)-3:find(round(K)==S0)+3];

% Truncate the outputed Greeks and strikes
K = K(ATM);
GreeksFFT = GreeksFFT(ATM,:);
 
% Adjust the FFT Greeks vega1 and vanna for v0
GreeksFFT(:,6:7) = GreeksFFT(:,6:7) .* (2*sqrt(v0));

%% Closed form Greeks and errors
GreeksCM = zeros(length(ATM),8);
[x w] = GenerateGaussLaguerre(32);

for k=1:length(K);
    for j=1:8
        GreeksCM(k,j)  = CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,tau,K(k),S0,r,v0,trap,Greek(j),x,w);
    end
end
for j=1:8
    error(j) = sum(abs(GreeksCM(:,j) - GreeksFFT(:,j)));
end

%% Write the results
fmtstring = repmat('%10.3f %8.3f',1,4);
fprintf('------------------------------------------------------------------------------------\n')
fprintf('                Call               Delta               Gamma              Theta     \n');
fprintf('Strike       FFT     CM         FFT     CM          FFT     CM         FFT     CM   \n');     
fprintf('------------------------------------------------------------------------------------\n')
for k=1:length(K)
    fprintf(['%6.2f ' fmtstring '\n'],K(k),...
            GreeksFFT(k,1),GreeksCM(k,1),...
            GreeksFFT(k,2),GreeksCM(k,2),...
            GreeksFFT(k,3),GreeksCM(k,3),...
            GreeksFFT(k,4),GreeksCM(k,4));
end
fprintf('Error   %16.3e %18.3e %18.3e %18.3e\n',error(1:4))
fmtstring = repmat('%10.3f %8.3f',1,4);
fprintf('------------------------------------------------------------------------------------\n')
fprintf('                 Rho               Vega1               Vanna              Volga     \n');
fprintf('Strike       FFT     CM         FFT     CM          FFT     CM         FFT     CM   \n');     
fprintf('------------------------------------------------------------------------------------\n')
for k=1:length(K)
    fprintf(['%6.2f ' fmtstring '\n'],K(k),...
            GreeksFFT(k,5),GreeksCM(k,5),...
            GreeksFFT(k,6),GreeksCM(k,6),...
            GreeksFFT(k,7),GreeksCM(k,7),...
            GreeksFFT(k,8),GreeksCM(k,8));
end
fprintf('Error   %16.3e %18.3e %18.3e %18.3e\n',error(1:4))
fprintf('------------------------------------------------------------------------------------\n')
