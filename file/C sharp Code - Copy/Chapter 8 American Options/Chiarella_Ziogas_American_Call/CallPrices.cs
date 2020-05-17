using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Chiarella_Ziogas_American_Call
{
    class CZPrices
    {
        public double CZEuroCall(double S0,double tau,HParam param,double K,double rf,double q,double[] Xs,double[] Ws)
        {
            Complex i = new Complex(0.0,1.0);
            double xi  = 0;
            Complex psi = 0;
            int Nx = Xs.Length;
            double[] Int1 = new double[Nx];
            double[] Int2 = new double[Nx];
            Complex I1,I2;

            // Compute the integrands
            CharFun CF = new CharFun();
            for(int k=0;k<=Nx-1;k++)
            {
                Complex phi = Xs[k];
                I1 = Complex.Exp(-i*phi*Math.Log(K)) * CF.CZCharFun(S0,tau,xi,param,K,rf,q,phi,psi,1) / (i*phi);
                Int1[k] = Ws[k] * I1.Real;
                I2 = Complex.Exp(-i*phi*Math.Log(K)) * CF.CZCharFun(S0,tau,xi,param,K,rf,q,phi,psi,2) / (i*phi);
                Int2[k] = Ws[k] * I2.Real;
            }

            // Define the probabilities
            double pi = Math.PI;
            double P1 = 0.5 + (1.0/pi)*Int1.Sum();
            double P2 = 0.5 + (1.0/pi)*Int2.Sum();

            // Return the call price
            return S0*Math.Exp(-q*tau)*P1 - K*Math.Exp(-rf*tau)*P2;
        }

        public double[] CZAmerCall(double S0,double tau,HParam param,double K,double rf,double q,double[] xs,double[] ws,double[] xt,double[] wt,
            int Nt,double b0,double b1,double a,double b,double c,double d,string DoubleType)
        {
            EarlyExercise EE = new EarlyExercise();
            double Euro    = CZEuroCall(S0,tau,param,K,rf,q,xs,ws);
            double Premium = EE.CZEarlyExercise(S0,tau,param,K,rf,q,xt,wt,xt,wt,Nt,b0,b1,a,b,c,d,DoubleType);
            double Amer    = Euro + Premium;
            double[] output = new double[2];
            output[0] = Amer;
            output[1] = Euro;
            return output;
        }
    }
}

