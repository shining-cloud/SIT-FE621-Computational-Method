function y = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,Greek)

% Double Heston (2009) call or put price by Trapezoidal Rule
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T = Time to maturity.
%   rf = Risk free rate.
%   q  = Dividend yield
%   param = Two sets of Double Heston parameters
%           [kappa1 theta1 sigma1 v01 rho1, 
%            kappa2 theta2 sigma2 v02 rho2] 
%   trap:  1 = "Little Trap" formulation (Gauthier and Possamai
%          0 = Original Heston formulation (Christoffersen et al)
%   b = Upper limit for Newton-Cotes
%   a = Lower limit for Newton-Cotes
%   N = number of integration points
% OUTPUT -------------------------------------------------------
%   The Double Heston call or put price

N = length(x);

for k=1:N;
	u = x(k);
	f2 = DoubleHestonCF(u  ,param,T,S,rf,q,trap);
	f1 = DoubleHestonCF(u-i,param,T,S,rf,q,trap);
    if strcmp(Greek,'Price')
        int2(k) = w(k) * real(exp(-i*u*log(K))/i/u*f2);
        int1(k) = w(k) * real(exp(-i*u*log(K))/i/u/exp((rf-q)*T)*f1/S);
    elseif strcmp(Greek,'Delta')
        int1(k) = w(k) * real(exp(-i*u*log(K))/i/u/exp((rf-q)*T)*f1/S);
    elseif strcmp(Greek,'Gamma')
        int1(k) = w(k) * real(exp(-i*u*log(K))/exp((rf-q)*T)*f1/S^2);
    elseif strcmp(Greek,'Rho')
        int2(k) = w(k) * real(exp(-i*u*log(K))/i/u*f2);
    elseif strcmp(Greek,'Vega11') || strcmp(Greek,'Vega12')
        if strcmp(Greek,'Vega11')
            v0choice = 1;
        else
            v0choice = 2;
        end
        B1 = B(u-i,param,T,trap,v0choice);
        B2 = B(u  ,param,T,trap,v0choice);
        df1 = f1*B1;
        df2 = f2*B2;
        int2(k) = w(k) * real(exp(-i*u*log(K))/i/u*df2);
        int1(k) = w(k) * real(exp(-i*u*log(K))/i/u/exp((rf-q)*T)*df1/S);
    elseif strcmp(Greek,'Volga1') || strcmp(Greek,'Volga2') 
        if strcmp(Greek,'Volga1')
            v0choice = 1;
            v0 = param(4);
        else
            v0choice = 2;
            v0 = param(9);
        end
        B1 = B(u-i,param,T,trap,v0choice);
        B2 = B(u  ,param,T,trap,v0choice);
        df1 = 2*B1*f1*(2*B1*v0 + 1);
        df2 = 2*B2*f2*(2*B2*v0 + 1);
        dint1(k) = w(k) * real(exp(-i*u*log(K))/i/u/exp((rf-q)*T)*df1/S);
        dint2(k) = w(k) * real(exp(-i*u*log(K))/i/u*df2);
    elseif strcmp(Greek,'Vanna1') || strcmp(Greek,'Vanna2')
        if strcmp(Greek,'Vanna1')
            v0choice = 1;
        else
            v0choice = 2;
        end
        B1 = B(u-i,param,T,trap,v0choice);
        df1 = f1*B1;
        dint1(k) = w(k) * real(exp(-i*u*log(K))/i/u/exp((rf-q)*T)*df1/S);
    elseif strcmp(Greek,'Theta')
        int2(k) = w(k) * real(exp(-i*u*log(K))/i/u*f2);
        int1(k) = w(k) * real(exp(-i*u*log(K))/i/u/exp((rf-q)*T)*f1/S);
        v01 = param(4);
        v02 = param(9);
        [dA1 dB11 dB21] = DiffTau(u-i,param,T,rf,q);
        df1 = f1*(dA1 + dB11*v01 + dB21*v02) - (rf-q)*f1;
        [dA2 dB12 dB22] = DiffTau(u  ,param,T,rf,q);
        df2 = f2*(dA2 + dB12*v01 + dB22*v02);
        dint2(k) = w(k) * real(exp(-i*u*log(K))/i/u*df2);
        dint1(k) = w(k) * real(exp(-i*u*log(K))/i/u/exp((rf-q)*T)*df1/S);
    end
end

% Output the price or the Greek
if strcmp(Greek,'Price')
    P1 = 1/2 + 1/pi*sum(int1);
    P2 = 1/2 + 1/pi*sum(int2);
	y = S*exp(-q*T)*P1 - K*exp(-rf*T)*P2;
elseif strcmp(Greek,'Delta')
    P1 = 1/2 + 1/pi*sum(int1);
    y = exp(-q*T)*P1;
elseif strcmp(Greek,'Gamma')
    y = exp(-q*T)*sum(int1)/pi;
elseif strcmp(Greek,'Rho')
    P2 = 1/2 + 1/pi*sum(int2);
    y = K*T*exp(-rf*T)*P2;
elseif strcmp(Greek,'Vega11') || strcmp(Greek,'Vega12')
    if strcmp(Greek,'Vega11')
        v0 = param(4);
    else
        v0 = param(9);
    end
    dP1 = 1/pi*sum(int1);
    dP2 = 1/pi*sum(int2);
	dC = S*exp(-q*T)*dP1 - K*exp(-rf*T)*dP2;
    y = dC*2*sqrt(v0);
elseif strcmp(Greek,'Volga1') || strcmp(Greek,'Volga2')
    dP1 = 1/pi*sum(dint1);
    dP2 = 1/pi*sum(dint2);
	y = S*exp(-q*T)*dP1 - K*exp(-rf*T)*dP2;
elseif strcmp(Greek,'Vanna1') || strcmp(Greek,'Vanna2')
    if strcmp(Greek,'Vanna1')
        v0 = param(4);
    else
        v0 = param(9);
    end
    dP1 = 1/pi*sum(dint1);
	y = 2*exp(-q*T)*sqrt(v0)*dP1;
elseif strcmp(Greek,'Theta')
    P1 = 1/2 + 1/pi*sum(int1);
    P2 = 1/2 + 1/pi*sum(int2);
    dP1 = 1/pi*sum(dint1);
    dP2 = 1/pi*sum(dint2);
    y = -S*exp(-q*T)*(-q*P1+dP1) + K*exp(-rf*T)*(-rf*P2+dP2);
end
