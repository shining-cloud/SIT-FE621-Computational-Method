function [bR bL LowerAP UpperAP LowerCF UpperCF] = FindLeeBounds(S,r,div,tau,param,trap,LoLimit,HiLimit)

% Finds the Roger Lee bounds on vol slope,
% and the Andersen and Piterbarg bounds on moments
% bR = Right bound on slope
% bL = Left bound on slope
% LowerAP, UpperAP = moment bounds using Andersen and Piterbarg
% LowerCF, UpperCF = moment bounds using the characteristic function

% Heston parameters
kappa = param(1);
theta = param(2);
sigma = param(3);
v0    = param(4);
rho   = param(5);

% Moment explosion setting
lambda = 1;

%% Plot to identify visually where the moments explode
%W =[LoLimit:.1:HiLimit];
% for j=1:length(W);
%     w = W(j);
%     T(j) = MomentExplode(w,lambda,sigma,kappa,rho);
% end
% plot(W,T)
% clear W T

%% Find the upper and lower moment bounds using Andersen and Piterbarg formula
% Upper moment bound
W = [1:HiLimit];
for k=1:15
    j = 1;
    T = Inf;
    while(T == inf)
        j = j+1;
        T = MomentExplode(W(j),lambda,sigma,kappa,rho);
    end
    e = 1/10^k;
    W = [W(j-1):e:W(j+1)];
    clear T
end
UpperAP = W(j);

% Lower moment bound
W = [LoLimit:0];
for k=1:15
    j = 1;
    T(j) = 0;
    while(T(j) < inf)
        j = j+1;
        T(j) = MomentExplode(W(j),lambda,sigma,kappa,rho);
    end
    e = 1/10^k;
    W = [W(j-1):e:W(j+1)];
    clear T
end
LowerAP = W(j);


%% Find moments using the CF
% Upper Moment.  Loop through until imag(CF) is encountered
CF = 0;
e = 1e-5;
W = 0.9*UpperAP;
while isreal(CF)
    phi = -i*W;
    CF = HestonCF(phi,kappa,theta,0,rho,sigma,tau,S,r,div,v0,trap);
    W = W+e;
end
UpperCF = W;

% Lower Moment.  Loop through until imag(CF) is encountered
CF = 0;
e = -1e-5;
W = 0.9*LowerAP;
while isreal(CF)
    phi = -i*W;
    CF = HestonCF(phi,kappa,theta,0,rho,sigma,tau,S,r,div,v0,trap);
    W = W+e;
end
LowerCF = W;

%% Roger Lee Moment formulas
p =  UpperAP - 1;
q = -LowerAP;

% Moment formula part 1 -- Limit of right tail slope
bR = 2 - 4*(sqrt(p^2 + p) - p);

% Moment formula part 2 -- Limit of the left tail slope
bL = 2 - 4*(sqrt(q^2 + q) - q);

