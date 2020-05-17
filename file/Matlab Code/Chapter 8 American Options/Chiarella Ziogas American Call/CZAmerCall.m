function [Amer Euro] = CZAmerCall(S0,tau,params,K,rf,q,xs,ws,xt,wt,Nt,b0,b1,a,b,c,d,DoubleType)

Euro    = CZEuroCall(S0,tau,params,K,rf,q,xs,ws);
Premium = CZEarlyExercise(S0,tau,params,K,rf,q,xt,wt,xt,wt,Nt,b0,b1,a,b,c,d,DoubleType);
Amer    = Euro + Premium;

