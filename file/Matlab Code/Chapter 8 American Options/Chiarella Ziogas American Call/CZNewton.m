function b = CZNewton(start,v0,v1,mat,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,B0,B1,gNum,tol,A,B,C,D,DoubleType)

% Increment for numerical derivative
db = 0.001;
diff = 1.1*tol;
b = start;
params0 = [kappa theta sigma v0 rho lambda];
params1 = [kappa theta sigma v1 rho lambda];
B = mat;

% Newton's method
while (abs(diff) > tol)
    if gNum == 1
        % First set of functions and derivative
        [CallA1 Euro] = CZAmerCall(exp(B0 + b*v0),mat,params0,K,r,q,xs,ws,xt,wt,Nt,B0,b,A,B,C,D,DoubleType);
        g0 = (log(CallA1 + K) - B0)/v0 - b;
        [CallA1 Euro] = CZAmerCall(exp(B0 + (b+db)*v0),mat,params0,K,r,q,xs,ws,xt,wt,Nt,B0,b+db,A,B,C,D,DoubleType);
        g  = (log(CallA1 + K) - B0)/v0 - (b+db);
        [CallA1 Euro] = CZAmerCall(exp(B0 + (b-db)*v0),mat,params0,K,r,q,xs,ws,xt,wt,Nt,B0,b-db,A,B,C,D,DoubleType);
        g_ = (log(CallA1 + K) - B0)/v0 - (b-db);
    else
        % Second set of functions and derivatives
        [CallA0 Euro] = CZAmerCall(exp(b + B1*v1),mat,params1,K,r,q,xs,ws,xt,wt,Nt,b,B1,A,B,C,D,DoubleType);
        g0 = (log(CallA0 + K) - v1*B1) - b;
        [CallA0 Euro] = CZAmerCall(exp(b+db + B1*v1),mat,params1,K,r,q,xs,ws,xt,wt,Nt,b+db,B1,A,B,C,D,DoubleType);
        g  = (log(CallA0 + K) - v1*B1) - (b+db);
        [CallA0 Euro] = CZAmerCall(exp(b-db + B1*v1),mat,params1,K,r,q,xs,ws,xt,wt,Nt,b-db,B1,A,B,C,D,DoubleType);
        g_ = (log(CallA0 + K) - v1*B1) - (b-db);
    end
    % The derivative
    dg = (g - g_)/2/db;
    % Newton's method
    b_new = b - g0/dg;
    diff = b_new - b;
    b = b_new;
end

