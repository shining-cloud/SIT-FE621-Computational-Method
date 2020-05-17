function [EuroPrice AmerPrice Euro Amer Yt V X Prob Branch] = BuildBivariateTree(S0,PutCall,Strike,T,rf,NT,kappa,theta,sigma,V0,rho,threshold)

% Constructs the bivariate Beliaeva-Nawalkha tree
% Uses cells for the node branches 
% Uses cells for the probabilities

dt = T/NT;
sigmay0 = sqrt(1-rho^2)*sqrt(V0);

%% Generate the volatility trees
[X V RBound M] = BuildVolTree(kappa,theta,sigma,V0,dt,NT,threshold);
if RBound ~= 2*NT-1
    disp('Warning: V-matrix not triangular')
end

% Find the column where the tree changes from triangular to rectangular
ColChange = RBound - M + 2;

%% The values of k(t) from Equation (11) and I(t) from (14)
k = zeros(2*NT-1,NT);
for t=1:NT
    for n=M-(t-1):M+(t-1);
        if V(n,t) > 0
            k(n,t) = ceil(sqrt(V(n,t)/V0));
        else
            k(n,t) = 1;
        end
    end
end


%% Find dimensions of the tree
maxK = max(k);
numY = [1 3];
numV = [1 3];
numRows = [1 9];
for t=3:NT
    numY(t) = numY(t-1) + 2*maxK(t);
    numV(t) = 2*t - 1;
    numRows(t) = numV(t)*numY(t);
end
NR = max(numRows);

%% Find the branch indices
Branch = cell(numRows(NT-1),NT-1);
Branch(1,1) = {[1:9]};
B = zeros(numV(NT-1)*numY(NT-1),9);

% To Branches
for t=2:NT-1
    nY = numY(t);
    First = maxK(t)+1;
    K = k(M-(t-1):M+(t-1),t);
    for n = 1:numRows(t);
        a = ceil(n/nY);
        b = mod(n-1,nY);
        % Find the middle-to-middle branches
        B(n,2) = First + (a-1)*numY(t+1) + b;
        B(n,5) = B(n,2) +      numY(t+1);
        B(n,8) = B(n,2) +    2*numY(t+1);
        % Find the rest of the branches
        B(n,1) = B(n,2) - K(a);
        B(n,3) = B(n,2) + K(a);
        B(n,4) = B(n,5) - K(a);
        B(n,6) = B(n,5) + K(a);
        B(n,7) = B(n,8) - K(a);
        B(n,9) = B(n,8) + K(a);
        Branch(n,t) = {B(n,:)};
    end
end
clear B;

% Adjust the last branch upward
for t=ColChange:NT-1
    for j=7:9
        Branch{numRows(t),t}(j) = Branch{numRows(t)-1,t}(j);
    end
end


%% Find the values for the indices and for the probabilities

% Log stock tree (Yt), stock price tree (St), and probabilities (Prob)
Yt = zeros(NR,NT);
Y0 = log(S0) - rho*V0/sigma;
Yt(1,1) = Y0;
Prob = cell(numRows(NT-1),NT-1);
X0 = X(M,1);

for t=1:NT-1;
    n = 0;
    J = -(t-1):(t-1);
    for j=1:numV(t);
        for r=1:numY(t);
            n = n + 1;
            Vt = V(M+J(j),t);
            Xt = X(M+J(j),t);
            Kt = k(M+J(j),t);
            NewBranch = cell2mat(Branch(n,t));
            muy = (rho*kappa/sigma - 0.5)*Vt;
            I = round(muy/Kt/sigmay0*sqrt(dt));
            if Yt(n,t) > 0;
                for s=1:9
                    if s==2 || s==5 || s==8
                        % Middle node
                        Yt(NewBranch(s),t+1) = Yt(n,t) + (I+0)*Kt*sigmay0*sqrt(dt);
                    elseif s==1 || s==4 || s==7
                        % Up node
                        Yt(NewBranch(s),t+1) = Yt(n,t) + (I+1)*Kt*sigmay0*sqrt(dt);
                    elseif s==3 || s==6 || s==9
                        % Down node
                        Yt(NewBranch(s),t+1) = Yt(n,t) + (I-1)*Kt*sigmay0*sqrt(dt);
                    end
                end
            else
                Branch(n,t) = {['No branch originates from (' num2str(t) ',' num2str(n) ')']};
            end
            if Yt(n,t) == 0;
                Prob(n,t) = {zeros(1,9)};
            else
                [pvu pvm pvd] = probV(Xt,X0,dt,kappa,theta,sigma);
                [pyu pym pyd] = probY(Vt,V0,Yt(n,t),dt,rho,sigma,kappa);
                prob = [pvu*pyu pvu*pym pvu*pyd pvm*pyu pvm*pym pvm*pyd pvd*pyu pvd*pym pvd*pyd];
                Prob(n,t) = {prob};
            end
        end
    end
end


%% Find the American and European option prices
Euro = zeros(NR,NT);
Amer = zeros(NR,NT);

% Last column for the stock price
t = NT;
ht = (rf - rho*kappa*theta/sigma)*(t-1)*dt;
n = 0;
J = -(t-1):(t-1);
ST = zeros(numRows(t),1);
for j=1:numV(t);
    for r=1:numY(t);
        n = n + 1;
        Vt = V(M+J(j),t);
        if Yt(n,t) > 0
            ST(n) = exp(Yt(n,t) + rho/sigma*Vt + ht);
        end
    end
end

% Payoff at maturity
if strcmp(PutCall,'C')
    Euro(:,NT) = max(ST - Strike, 0);
    Amer(:,NT) = max(ST - Strike, 0);
elseif strcmp(PutCall,'P')
    Euro(:,NT) = max(Strike - ST, 0);
    Amer(:,NT) = max(Strike - ST, 0);
end

for t = NT-1:-1:1
    n = 0;
    ht = (rf - rho*kappa*theta/sigma)*(t-1)*dt;
    J = -(t-1):(t-1);
    for j=1:numV(t);
        for r=1:numY(t);
            n = n + 1;
            Vt = V(M+J(j),t);
            if Yt(n,t) > 0
                St = exp(Yt(n,t) + rho/sigma*Vt + ht);
                P = cell2mat(Prob(n,t));
                B = cell2mat(Branch(n,t));
                Euro(n,t) = P(1)*Euro(B(1),t+1) + P(2)*Euro(B(2),t+1) + P(3)*Euro(B(3),t+1)...
                    + P(4)*Euro(B(4),t+1) + P(5)*Euro(B(5),t+1) + P(6)*Euro(B(6),t+1)...
                    + P(7)*Euro(B(7),t+1) + P(8)*Euro(B(8),t+1) + P(9)*Euro(B(9),t+1);
                Euro(n,t) = exp(-rf*dt)*Euro(n,t);
                Amer(n,t) = P(1)*Amer(B(1),t+1) + P(2)*Amer(B(2),t+1) + P(3)*Amer(B(3),t+1)...
                    + P(4)*Amer(B(4),t+1) + P(5)*Amer(B(5),t+1) + P(6)*Amer(B(6),t+1)...
                    + P(7)*Amer(B(7),t+1) + P(8)*Amer(B(8),t+1) + P(9)*Amer(B(9),t+1);
                Amer(n,t) = exp(-rf*dt)*Amer(n,t);
            end
            if strcmp(PutCall,'C')
                Amer(n,t) = max(St - Strike, Amer(n,t));
            elseif strcmp(PutCall,'P')
                Amer(n,t) = max(Strike - St, Amer(n,t));
            end
        end
    end
end

EuroPrice = Euro(1,1);
AmerPrice = Amer(1,1);

