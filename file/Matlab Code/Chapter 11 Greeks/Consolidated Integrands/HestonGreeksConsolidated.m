function y = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,Greek)

N = length(x);
for j=1:N
    phi = x(j);
    f1 = HestonCF(phi-i,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap);
    f2 = HestonCF(phi,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap);
    if strcmp(Greek,'Price')
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (f1 - K*f2));
    elseif strcmp(Greek,'Delta')
        df1 = f1 * (i*phi+1)/S;
        df2 = f2 * i*phi/S;
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (df1 - K*df2));
    elseif strcmp(Greek,'Gamma')
        df1 = f1 * (-phi^2+i*phi)/S^2;
        df2 = f2 * (-phi^2-i*phi)/S^2;
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (df1 - K*df2));
    elseif strcmp(Greek,'Theta')
        [dC1 dD1] = DiffTau(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        [dC2 dD2] = DiffTau(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        df1 = f1*(dC1+dD1*v0);
        df2 = f2*(dC2+dD2*v0);
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (-r*(f1-K*f2)+df1-K*df2));
    elseif strcmp(Greek,'Rho')
        df1 = f1 * i*(phi-i)*tau;
        df2 = f2 * i*phi*tau;
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (-tau*(f1-K*f2)+df1-K*df2));
    elseif strcmp(Greek,'Vega1')
        D1 = D(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        D2 = D(phi  ,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        df1 = f1*D1;
        df2 = f2*D2;
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (df1 - K*df2));
    elseif strcmp(Greek,'Vanna')
        D1 = D(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        D2 = D(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        df1 = f1 * D1*(i*phi+1)/S;
        df2 = f2 * D2*i*phi/S;
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (df1 - K*df2));
    elseif strcmp(Greek,'Volga')
        D1 = D(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        D2 = D(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
        df1 = f1*D1;
        df2 = f2*D2;
        d2f1 = f1*D1*D1;
        d2f2 = f2*D2*D2;
        int(j) = w(j) * real(exp(-i*phi*log(K)-r*tau)/i/phi * (4*v0*(d2f1-K*d2f2) + 2*(df1-K*df2)));
    end
end
Integral = sum(int);

if strcmp(Greek,'Price')
    y = S*exp(-q*tau)/2 - K*exp(-r*tau)/2 + (1/pi)*Integral;
elseif strcmp(Greek,'Delta')
    y = exp(-q*tau)/2 + (1/pi)*Integral;
elseif strcmp(Greek,'Theta')
    y = S/2*q*exp(-q*tau) - r/2*K*exp(-r*tau) - (1/pi)*Integral;
elseif strcmp(Greek,'Rho')
    y = tau/2*K*exp(-r*tau) + (1/pi)*Integral;
elseif strcmp(Greek,'Gamma') || strcmp(Greek,'Volga')
    y = (1/pi)*Integral;
elseif strcmp(Greek,'Vega1') || strcmp(Greek,'Vanna')
    y = (1/pi)*Integral*2*sqrt(v0);
end

