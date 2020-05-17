using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace MultiDomain_Integration_Zhu
{
    class NewtonCotesPrice
    {
        // Heston price by Newton Coates quadratures
        public double HestonPriceNewtonCotes(HParam param,OpSet settings,int method,double a,double b,int N)
        {
            // Built the integration grid
            // For Simpson's 3/8 rule, the code ensures that N-1 is divisible by 3
            double h = (b-a)/(Convert.ToDouble(N)-1.0);
            double[] phi = new double[N];
            phi[0] = a;
            phi[N-1] = b;
            for(int k=1;k<=N-2;k++)
                phi[k] = phi[k-1] + h;

            double[] int1 = new double[N];
            double[] int2 = new double[N];

            //   Integration methods
            HestonPrice HP = new HestonPrice();
            if(method==1)
            {
                // Mid-Point rule --------------------------------------
                double[] wt = new double[N];
                for(int k=0;k<=N-1;k++)
                    wt[k] = h;
                for(int k=0;k<=N-2;k++)
                {
                    double mid = (phi[k] + phi[k+1])/2.0;
                    int1[k] = wt[k] * HP.HestonProb(mid,param,settings,1);
                    int2[k] = wt[k] * HP.HestonProb(mid,param,settings,2);
                }
                int1[N-1] = 0.0;
                int2[N-1] = 0.0;
            }
            else if(method==2)
            {
                // Trapezoidal rule --------------------------------------------------
                double[] wt = new double[N];
                wt[0]   = 0.5*h;
                wt[N-1] = 0.5*h;
                for(int k=1;k<=N-2;k++)
                    wt[k] = h;
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = wt[k] * HP.HestonProb(phi[k],param,settings,1);
                    int2[k] = wt[k] * HP.HestonProb(phi[k],param,settings,2);
                }
            }
            else if(method==3)
            {
                // Simpson's Rule ----------------------------------------------------
                double[] wt = new double[N];
                wt[0]   = h/3.0;
                wt[N-1] = h/3.0;
                for(int k=1;k<=N-1;k++)
                    wt[k] = (h/3.0) * (3 + Math.Pow(-1,k+1));
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = wt[k] * HP.HestonProb(phi[k],param,settings,1);
                    int2[k] = wt[k] * HP.HestonProb(phi[k],param,settings,2);
                }

            }
            else if(method==4)
            {
                // Simpson's 3/8 rule --------------------------------------------
                // Ensure that N-1 is divisible by 3
                N = N - (N % 3);

                // Build the new grid
                h = (b-a)/(N-1.0);
                double[] phi2 = new double[N];
                phi2[0] = a;
                phi2[N-1] = b;
                for(int k=1;k<=N-2;k++)
                    phi2[k] = phi2[k-1] + h;

                double[] wt = new double[N];
                wt[0]   = 3.0/8.0*h;
                wt[1]   = 9.0/8.0*h;
                wt[2]   = 9.0/8.0*h;
                wt[N-1] = 3.0/8.0*h;
                for(int k=3;k<=N-1;k++)
                {
                    if((k % 3) == 1)
                        wt[k] = 6.0/8.0*h;
                    else
                        wt[k] = 9.0/8.0*h;
                }
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = wt[k] * HP.HestonProb(phi2[k],param,settings,1);
                    int2[k] = wt[k] * HP.HestonProb(phi2[k],param,settings,2);
                }
            }

            // The integrals
            double I1 = int1.Sum();
            double I2 = int2.Sum();

            // The probabilities P1 and P2
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*I1;
            double P2 = 0.5 + 1.0/pi*I2;

            // The call price
            double S = settings.S;
            double K = settings.K;
            double q = settings.q;
            double r = settings.r;
            double T = settings.T;
            string PutCall = settings.PutCall;

            double HestonC = S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);

            // Output the option price
            if(PutCall == "C")
                return HestonC;
            else
                return HestonP;
        }
    }
}
