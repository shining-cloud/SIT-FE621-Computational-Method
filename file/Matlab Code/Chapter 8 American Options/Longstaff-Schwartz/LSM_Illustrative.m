function [EuroPrice AmerPrice] = LSM_Illustrative(S,K,r,q,T,NT,NS,PutCall,XmatrixHandle)

% Longstaff-Schwartz for American puts and calls
% Also returns Euro price

dt = T/NT;              % Time increment.

% Initialize the Cash Flows.
CF = zeros(NS,NT);

% Set the last cash flows to the intrinsic value.
if strcmp(PutCall,'P')
    CF(:,NT) = max(K - S(:,NT), 0);
elseif strcmp(PutCall,'C')
    CF(:,NT) = max(S(:,NT) - K, 0);
end

% Number of columns in the X matrix
NX = length(XmatrixHandle);

% European price
EuroPrice = exp(-r*T)*mean(CF(:,NT));

% Work backwards through the stock prices until time t=2.
% We could work through to time t=1 but the regression will not be
% of full rank at time 1, so this is cleaner.
for t = NT-1:-1:2
    % Indices for stock paths in-the-money at time t
    if strcmp(PutCall,'P')
        I = find(S(:,t) < K);
    else
        I = find(S(:,t) > K);
    end
    fprintf('t = %1.0f --------------------------------------------------------------\n',t)
    
    % Stock paths in-the-money at time t
    X = S(I,t);
    fmtstr1 = ['%7.0f' repmat('%10.0f ',1,length(I)-1)];
    fprintf(['In the money indices         ' fmtstr1 '\n'],I);
    fmtstr2 = ['%7.4f' repmat('%10.4f ',1,length(I)-1)];
    fprintf(['In the money prices          ' fmtstr2 '\n'],X);

    % Cash flows at time t+1, discounted one period
    Y = CF(I,t+1)*exp(-r*dt);
%   fprintf(['Time t+1 cash flows          ' fmtstr2 '\n'],Y);

    % Design matrix for regression to predict cash flows
    Z = zeros(length(X),NX);
    for k=1:NX
        Z(:, k) = feval(XmatrixHandle{k}, X);
    end
    
    % Regression parameters and predicted cash flows
    beta = Z\Y;
    PredCF = Z*beta;
    fprintf(['Predicted cash flows         ' fmtstr2 '\n'],PredCF);

    % Indices for stock paths where immediate exercise is optimal
    if strcmp(PutCall,'P')
        J = max(K - X, 0) > PredCF;
    else
        J = max(X - K, 0) > PredCF;
    end

    % In-the-money stock price indices where immediate exercise is optimal
    fprintf(['Immediate Exercise           ' fmtstr2 '\n'],K-X');
    Ex = I(J);
    fmtstr3 = ['%7.0f' repmat('%10.0f ',1,length(Ex)-1)];
    fprintf(['Paths of optimal exercise    ' fmtstr3 '\n'],Ex);
   

    % All other stock path indices --> continuation is optimal
    Co = setdiff((1:NS),Ex)';

    % Replace cash flows with exercise value where exercise is optimal
    if strcmp(PutCall,'P')
        CF(Ex,t) = max(K - X(J), 0);
    else
        CF(Ex,t) = max(X(J) - K, 0);
    end
%     fmtstr4 = ['%7.4f' repmat('%10.4f ',1,length(Ex)-1)];
%     fprintf(['Updated cash flows           ' fmtstr4 '\n'],CF(Ex,t))

    fmtstr4 = ['%7.4f' repmat('%10.4f ',1,length(I)-1)];
    fprintf(['Updated cash flows           ' fmtstr4 '\n'],CF(I,t))
    
    % Continued CF are discounted back one period.
    CF(Co,t) = exp(-r*dt)*CF(Co,t+1);
end

% The American option price = cash flow at period 2 discounted one peroid
AmerPrice = exp(-r*dt)*mean(CF(:,2));

