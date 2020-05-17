function [EuroPrice AmerPrice Euro Amer Yt V X Prob Branch] = BuildBivariateTree3(S0,PutCall,Strike,T,rf,NT,kappa,theta,sigma,V0,rho,threshold)

% Constructs the bivariate Beliaeva-Nawalkha tree
% Uses matrices for the node branches 
% Uses matrices for the probabilities

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
Branch = zeros(numRows(NT-1),9*(NT-1));
for j=1:9
    Branch(1,j) = j;
end

for t=2:NT-1
    nY = numY(t);
    First = maxK(t)+1;
    K = k(M-(t-1):M+(t-1),t);
    for n = 1:numRows(t);
        a = ceil(n/nY);
        b = mod(n-1,nY);
        % Find the middle-to-middle branches
        Branch(n,9*(t-1)+2) = First + (a-1)*numY(t+1) + b;
        Branch(n,9*(t-1)+5) = Branch(n,9*(t-1)+2) +      numY(t+1);
        Branch(n,9*(t-1)+8) = Branch(n,9*(t-1)+2) +    2*numY(t+1);
        % Find the rest of the branches
        Branch(n,9*(t-1)+1) = Branch(n,9*(t-1)+2) - K(a);
        Branch(n,9*(t-1)+3) = Branch(n,9*(t-1)+2) + K(a);
        Branch(n,9*(t-1)+4) = Branch(n,9*(t-1)+5) - K(a);
        Branch(n,9*(t-1)+6) = Branch(n,9*(t-1)+5) + K(a);
        Branch(n,9*(t-1)+7) = Branch(n,9*(t-1)+8) - K(a);
        Branch(n,9*(t-1)+9) = Branch(n,9*(t-1)+8) + K(a);
    end
end

% Adjust the last branch upward
for t=ColChange:NT-1
    Branch(numRows(t),9*(t-1)+7:9*(t-1)+9) = Branch(numRows(t)-1,9*(t-1)+7:9*(t-1)+9);
end


%% Find the values for the indices and for the probabilities
 
% Log stock tree (Yt), stock price tree (St), and probabilities (Prob)
Yt = zeros(NR,NT);
Y0 = log(S0) - rho*V0/sigma;
Yt(1,1) = Y0;
Prob = zeros(numRows(NT-1),9*(NT-1));
X0 = X(M,1);

for t=1:NT-1;
    n = 0;
    J = -(t-1):(t-1);
    for j=1:numV(t);
        for r=1:numY(t)
            n = n + 1;
            Vt = V(M+J(j),t);
            Xt = X(M+J(j),t);
            Kt = k(M+J(j),t);
            NewBranch = Branch(n,9*(t-1)+1:9*(t-1)+9);
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
            end
            if Yt(n,t) > 0
                [pvu pvm pvd] = probV(Xt,X0,dt,kappa,theta,sigma);
                [pyu pym pyd] = probY(Vt,V0,Yt(n,t),dt,rho,sigma,kappa);
                prob = [pvu*pyu pvu*pym pvu*pyd pvm*pyu pvm*pym pvm*pyd pvd*pyu pvd*pym pvd*pyd];
                Prob(n,9*(t-1)+1:9*(t-1)+9) = prob;
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
                P = Prob(n,9*(t-1)+1:9*(t-1)+9);
                B = Branch(n,9*(t-1)+1:9*(t-1)+9);
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

