using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Benhamou_Gobet_Miri_Piecewise_Parameters
{
    class HestonPriceTD
    {
        // Mikhailov and Nogel (2003) Time dependent Heston Integrand
        public double HestonProbTD(double phi,HParam param,double[,] param0,double tau,double[] tau0,OpSet settings,double K,int Pnum)
        {
            CDcoefficients CD = new CDcoefficients();

            Complex i = new Complex(0.0,1.0);                   // Imaginary unit
            double S = settings.S;
            double r = settings.r;
            double q = settings.q;
            double x  = Math.Log(S);
            Complex C = new Complex(0.0,0.0);
            Complex D = new Complex(0.0,0.0);
            Complex C0 = new Complex(0.0,0.0);
            Complex D0 = new Complex(0.0,0.0);

            int N = tau0.Length;
            double v0,T;

            // Create the past Cj and Dj values corresponding to the old maturities
            for(int t=0;t<=N-1;t++)
            {
                HParam param2 = new HParam();
                param2.kappa = param0[t,0];
                param2.theta = param0[t,1];
                param2.sigma = param0[t,2];
                param2.v0    = param0[t,3];
                param2.rho   = param0[t,4];
                T = tau0[t];
                if(t==0)
                {
                    D0 = 0.0;
                    C0 = 0.0;
                }
                else
                {
                    D0 = D;
                    C0 = C;
                }
                C = CD.Ct(phi,param2,r,q,T,settings,Pnum,C0,D0);
                D = CD.Dt(phi,param2,r,q,T,settings,Pnum,C0,D0);
            }

            // Current Cj and Dj values
            v0 = param.v0;
            T = tau;
            if(tau0[0] == 0.0)
            {
                D0 = 0.0;
                C0 = 0.0;
            }
            else
            {
                D0 = D;
                C0 = C;
            }
            C = CD.Ct(phi,param,r,q,T,settings,Pnum,C0,D0);
            D = CD.Dt(phi,param,r,q,T,settings,Pnum,C0,D0);

            // The characteristic function.
            Complex f = Complex.Exp(C + D*v0 + i*phi*x);

            // The integrand.
            Complex integrand = Complex.Exp(-i*phi*Complex.Log(K))*f/i/phi;

            // Return the real part of the integrand
            return integrand.Real;
        }

        // Mikhailov and Nogel (2003) Heston Price by Gauss-Laguerre Integration
        public double HestonPriceGaussLaguerreTD(HParam param, double[,] param0,double tau,double[] tau0,OpSet settings,double K,double[] x, double[] w)
        {
            double[] int1 = new Double[32];
            double[] int2 = new Double[32];

            // Numerical integration
            for (int j=0; j<=31; j++)
            {
                int1[j] = w[j] * HestonProbTD(x[j],param,param0,tau,tau0,settings,K,1);
                int2[j] = w[j] * HestonProbTD(x[j],param,param0,tau,tau0,settings,K,2);
            }

            // Define P1 and P2
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*int1.Sum();
            double P2 = 0.5 + 1.0/pi*int2.Sum();

            // The call price
            double S = settings.S;
            double r = settings.r;
            double q = settings.q;
            string PutCall = settings.PutCall;
            double HestonC = S*Math.Exp(-q*tau)*P1 - K*Math.Exp(-r*tau)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*tau) + K*Math.Exp(-r*tau);

            // Output the option price
            if (PutCall == "C")
	            return HestonC;
            else
	            return HestonP;
        }
    }
}

