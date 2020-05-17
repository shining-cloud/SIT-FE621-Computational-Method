using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Chiarella_Ziogas_American_Call
{
    class DoubleIntegral
    {
        // Double integral using Gauss-Legendre
        public double DoubleGaussLegendre(double S0,double tau,HParam param,double K,double rf,double q,double b0,double b1,
            double[] xt,double[] wt,double[] xs,double[] ws,double a,double b,double c,double d,int funNum)
        {
            Complex i = new Complex(0.0,1.0);
            int Nt = xt.Length;
            int Ns = xs.Length;
            double h1 = (b-a)/2.0;
            double h2 = (b+a)/2.0;
            double k1 = (d-c)/2.0;
            double k2 = (d+c)/2.0;
            double time,phi,realfun,qr=0.0;
            Complex fun;

            // Select rate or dividend
            if(funNum == 1)
                qr = q;
            else if(funNum == 2)
                qr = rf;
            double Int = 0.0;

            // Double integral
            CharFun CF = new CharFun();
            for(int t=0;t<=Nt-1;t++)
            {
                time = h1*xt[t] + h2;
                for(int x=0;x<=Ns-1;x++)
                {
                    phi = k1*xs[x] + k2;
                    fun = Complex.Exp(-b0*i*phi) * CF.CZCharFun(S0,tau,time,param,K,rf,q,phi,-b1*phi,funNum)/(i*phi);
                    realfun = Math.Exp(qr*time) * fun.Real;
                    Int += h1*k1*wt[t]*ws[x] * realfun;
                }
            }
            return Int;
        }

        // Double integral using trapezoidal rule
        public double DoubleTrapezoidal(HParam param,double S0,double K,double tau,double rf,double q,
            double b0,double b1,double[] X,double[] T,int funNum)
        {
            Complex i = new Complex(0.0,1.0);
            double a,b,c,d,qr=0.0;
            Complex fun1,fun2,fun3,fun4;
            double g1,g2,g3,g4,term1,term2,term3;
            int Nt = T.Length;
            int Nx = X.Length;

            // Select rate or dividend
            if(funNum == 1)
                qr = q;
            else if(funNum == 2)
                qr = rf;
            double[,] Int = new double[Nt,Nx];
            double sumInt = 0.0;

            // Double trapezoidal rule
            CharFun CF = new CharFun();
            for(int t=1;t<=Nt-1;t++)
            {
                a = T[t-1];
                b = T[t];
                for(int x=1;x<=Nx-1;x++)
                {
                    c = X[x-1];
                    d = X[x];
                    fun1 = Complex.Exp(-b0*i*c) * CF.CZCharFun(S0,tau,a,param,K,rf,q,c,-b1*c,funNum) / (i*c);
                    fun2 = Complex.Exp(-b0*i*d) * CF.CZCharFun(S0,tau,a,param,K,rf,q,d,-b1*d,funNum) / (i*d);
                    fun3 = Complex.Exp(-b0*i*c) * CF.CZCharFun(S0,tau,b,param,K,rf,q,c,-b1*c,funNum) / (i*c);
                    fun4 = Complex.Exp(-b0*i*d) * CF.CZCharFun(S0,tau,b,param,K,rf,q,d,-b1*d,funNum) / (i*d);
                    g1 = Math.Exp(qr*a) * fun1.Real;
                    g2 = Math.Exp(qr*a) * fun2.Real;
                    g3 = Math.Exp(qr*b) * fun3.Real;
                    g4 = Math.Exp(qr*b) * fun4.Real;
                    term1 = g1 + g2 + g3 + g4;
                    fun1 = Complex.Exp(-b0*i*c) * CF.CZCharFun(S0,tau,(a+b)/2.0,param,K,rf,q,c,-b1*c,funNum) / (i*c);
                    fun2 = Complex.Exp(-b0*i*d) * CF.CZCharFun(S0,tau,(a+b)/2.0,param,K,rf,q,d,-b1*d,funNum) / (i*d);
                    fun3 = Complex.Exp(-b0*i*(c+d)/2.0) * CF.CZCharFun(S0,tau,a,param,K,rf,q,(c+d)/2.0,-b1*(c+d)/2.0,funNum) / (i*(c+d)/2.0);
                    fun4 = Complex.Exp(-b0*i*(c+d)/2.0) * CF.CZCharFun(S0,tau,b,param,K,rf,q,(c+d)/2.0,-b1*(c+d)/2.0,funNum) / (i*(c+d)/2.0);
                    g1 = Math.Exp(qr*(a+b)/2.0) * fun1.Real;
                    g2 = Math.Exp(qr*(a+b)/2.0) * fun2.Real;
                    g3 = Math.Exp(qr*a) * fun3.Real;
                    g4 = Math.Exp(qr*b) * fun4.Real;
                    term2 = g1 + g2 + g3 + g4;
                    fun1 = Complex.Exp(-b0*i*(c+d)/2.0) * CF.CZCharFun(S0,tau,(a+b)/2.0,param,K,rf,q,(c+d)/2.0,-b1*(c+d)/2.0,funNum) / (i*(c+d)/2.0);
                    term3 = Math.Exp(qr*(a+b)/2.0) * fun1.Real;
                    Int[t,x] = (b-a)*(d-c)/16.0*(term1 + 2*term2 + 4*term3);
                    sumInt += Int[t,x];
                }
            }
            return sumInt;
        }
    }
}
