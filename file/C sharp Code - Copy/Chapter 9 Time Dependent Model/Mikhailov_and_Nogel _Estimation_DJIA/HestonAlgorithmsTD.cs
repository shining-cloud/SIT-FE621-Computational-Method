using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Mikhailov_and_Nogel_Estimation_DJIA
{
    class HestonPriceTD
    {
        // Time dependent Heston Integrand
        public double MNProb(double phi,HParam param,double[,] param0,double tau,double[] tau0,double S,double K,double r,double q,int Pnum,int Trap)
        {
            CDCoefficients CD = new CDCoefficients();

            Complex i = new Complex(0.0,1.0);                   // Imaginary unit
            double x  = Math.Log(S);
            Complex C = new Complex(0.0,0.0);
            Complex D = new Complex(0.0,0.0);
            Complex C0 = new Complex(0.0,0.0);
            Complex D0 = new Complex(0.0,0.0);

            int N = tau0.Length;
            double v0,T;

            // Create the past Cj and Dj values corresponding to the old maturities
            if(tau0[0] != 0.0)
            {
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
                    C = CD.Ct(phi,param2,S,K,r,q,T,Pnum,Trap,C0,D0);
                    D = CD.Dt(phi,param2,S,K,r,q,T,Pnum,Trap,C0,D0);
                }
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
            C = CD.Ct(phi,param,S,K,r,q,T,Pnum,Trap,C0,D0);
            D = CD.Dt(phi,param,S,K,r,q,T,Pnum,Trap,C0,D0);

            // The characteristic function.
            Complex f = Complex.Exp(C + D*v0 + i*phi*x);

            // The integrand.
            Complex integrand = Complex.Exp(-i*phi*Complex.Log(K))*f/i/phi;

            // Return the real part of the integrand
            return integrand.Real;
        }

        // Heston Price by Gauss-Laguerre Integration
        public double MNPriceGaussLaguerre(HParam param,double[,] param0,double tau,double[] tau0,double S,double K,double r,double q,string PutCall,int Trap,double[] x,double[] w)
        {
            double[] int1 = new Double[32];
            double[] int2 = new Double[32];

            // Numerical integration
            for (int j=0; j<=31; j++)
            {
                int1[j] = w[j] * MNProb(x[j],param,param0,tau,tau0,S,K,r,q,1,Trap);
                int2[j] = w[j] * MNProb(x[j],param,param0,tau,tau0,S,K,r,q,2,Trap);
            }

            // Define P1 and P2
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*int1.Sum();
            double P2 = 0.5 + 1.0/pi*int2.Sum();

            // The call price
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

