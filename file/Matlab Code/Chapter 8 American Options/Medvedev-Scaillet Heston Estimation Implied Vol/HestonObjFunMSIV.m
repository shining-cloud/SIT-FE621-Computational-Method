function y = HestonObjFunMSIV(param,S,rf,q,K,T,trap,MktIV,method,A,B,N,yinf,NumTerms,a,b,Tol,MaxIter)

% Heston parameters
kappa  = param(1); 
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);

[NK,NT] = size(MktIV);
ModelPrice = zeros(NK,NT);
ModelIV    = zeros(NK,NT);
error      = zeros(NK,NT);

for k=1:NK
	for t=1:NT
         ModelPrice(k,t) = MSPrice(S,K(k),T(t),rf,q,param,trap,method,A,B,N,NumTerms,yinf);
         ModelIV(k,t) = BisecMSIV(S,K(k),rf,q,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
         error(k,t) = (ModelIV(k,t) - MktIV(k,t))^2 / MktIV(k,t);
	end
end

% MSE loss function
y = sum(sum(error));

fprintf('%10.4e %10.4f %10.4f %10.4f %10.4f %10.4f \n',y,kappa,theta,sigma,v0,rho);


