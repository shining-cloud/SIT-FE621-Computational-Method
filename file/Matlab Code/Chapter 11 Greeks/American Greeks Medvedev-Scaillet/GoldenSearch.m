function y = GoldenSearch(a,b,tol,MaxIter,theta,K,params,r,q,T,NumTerms)

% Golden Ratio
GR = (sqrt(5) - 1)/2;
k = 0;

% Implement the Golden Section Search Method until tolerance is reached
tic
while (abs(b-a) > tol)
    k = k+1;
    x1 = a;
    x2 = a + (1-GR)*(b-a);
    x3 = a + GR*(b-a);
    x4 = b;
    f1 = -MSPutHeston(x1,theta,K,params,r,q,T,NumTerms);
    f2 = -MSPutHeston(x2,theta,K,params,r,q,T,NumTerms);
    f3 = -MSPutHeston(x3,theta,K,params,r,q,T,NumTerms);
    f4 = -MSPutHeston(x4,theta,K,params,r,q,T,NumTerms);
    if (f1 > f2) && (f2 < f3)
        b = x3;
    elseif (f2 > f3) && (f3 < f4);
        a = x2;
    end;
    if k > MaxIter;   % Break if k > MaxIterations
        break;
    end
end;
time = toc;

% Return the average of endpoints
y = (a+b)/2;
fprintf('Required %d iterations and %5.2f seconds \n',k,time);

