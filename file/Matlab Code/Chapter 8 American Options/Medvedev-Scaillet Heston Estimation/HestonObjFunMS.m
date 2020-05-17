function y = HestonObjFunMS(param,S,rf,q,K,T,trap,MktPrice,method,A,B,N,lb,ub,yinf,NumTerms)

kappa  = param(1); 
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);

[NK,NT] = size(MktPrice);
ModelPrice = zeros(NK,NT);
error      = zeros(NK,NT);

for k=1:NK
	for t=1:NT
         ModelPrice(k,t) = MSPrice(S,K(k),T(t),rf,q,param,trap,method,A,B,N,NumTerms,yinf);
         error(k,t) = (ModelPrice(k,t) - MktPrice(k,t))^2;
	end
end
% MSE loss function
y = sum(sum(error)) / (NT*NK);

% Impose penalty for inadmissible parameter values
if kappa<=lb(1)   || theta<=lb(2) || sigma<=lb(3) || v0<=lb(4) || rho<=lb(5) || ...
   kappa>=ub(1)   || theta>=ub(2) || sigma>=ub(3) || v0>=ub(4) || rho>=ub(5)
    y = 1e100;
end

fprintf('%10.4e %10.4f %10.4f %10.4f %10.4f %10.4f \n',y,kappa,theta,sigma,v0,rho);

