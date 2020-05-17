using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;


namespace Chiarella_Ziogas_American_Call
{
    class EarlyExercise
    {
        // The early exercise price using the composite trapezoidal rule

        public double CZEarlyExercise(double S0,double tau,HParam param,double K,double rf,double q,
            double[] xs,double[] ws,double[] xt,double[] wt,int Nt,double b0,double b1,
            double a,double b,double c,double d,string DoubleType)
        {
            DoubleIntegral DI = new DoubleIntegral();
            double Int1=0.0, Int2=0.0;
            if(DoubleType == "GLe")
            {
                Int1 = DI.DoubleGaussLegendre(S0,tau,param,K,rf,q,b0,b1,xt,wt,xt,wt,a,b,c,d,1);
                Int2 = DI.DoubleGaussLegendre(S0,tau,param,K,rf,q,b0,b1,xt,wt,xt,wt,a,b,c,d,2);
            }
            else if(DoubleType == "Trapz")
            {
                double ht = (b-a)/Nt;
                double hs = (d-c)/Nt;
                double[] T = new double[Nt+1];
                double[] X = new double[Nt+1];
                for(int j=0;j<=Nt;j++)
                {
                    T[j] = a + j*ht;
                    X[j] = c + j*hs;
                }

                Int1 = DI.DoubleTrapezoidal(param,S0,K,tau,rf,q,b0,b1,X,T,1);
                Int2 = DI.DoubleTrapezoidal(param,S0,K,tau,rf,q,b0,b1,X,T,2);
            }
            else
                return 0.0;

            // The early exercise premium
            double pi = Math.PI;
            double V1 = S0*(1-Math.Exp(-q*tau))/2 + (1/pi)*S0*q*Math.Exp(-q*tau)*Int1;
            double V2 = K*(1-Math.Exp(-rf*tau))/2 + (1/pi)*K*rf*Math.Exp(-rf*tau)*Int2;

            return V1 - V2;
        }
    }
}
