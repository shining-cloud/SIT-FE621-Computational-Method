using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Risk_Neutral_Density
{
    class HestonPriceMD
    {
        // Heston Price by Gauss-Legendre Integration
        public OutputMD HestonPriceGaussLegendreMD(HParam param,OpSet settings,double[] xGLe,double[] wGLe,double[] A,double tol)
        {
            HestonPrice HP = new HestonPrice();
            int nA = A.Length;
            int nX = xGLe.Length;
            double[,] int1 = new double[nA,nX];
            double[,] int2 = new double[nA,nX];
            double[] sum1 = new double[nA];
            double[] sum2 = new double[nA];

            // Numerical integration
            int nj = 0;
            for(int j=1;j<=nA-1;j++)
            {
                // Counter for the last point of A used
                nj += 1;
                sum1[j] = 0.0;
                for(int k=0;k<=nX-1;k++)
                {
                    // Lower and upper and limits of the subdomain
                    double a = A[j-1];
                    double b = A[j];
                    double X = (a+b)/2.0 + (b-a)/2.0*xGLe[k];
                    int1[j,k] = wGLe[k] * HP.HestonProb(X,param,settings,1)*(b-a)/2.0;
                    int2[j,k] = wGLe[k] * HP.HestonProb(X,param,settings,2)*(b-a)/2.0;
                    sum1[j] += int1[j,k];
                    sum2[j] += int2[j,k];
                }
                if(Math.Abs(sum1[j])<tol && Math.Abs(sum2[j])<tol)
                    break;
            }

            // Define P1 and P2
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*sum1.Sum();
            double P2 = 0.5 + 1.0/pi*sum2.Sum();

            // The call price
            double S = settings.S;
            double K = settings.K;
            double T = settings.T;
            double r = settings.r;
            double q = settings.q;
            string PutCall = settings.PutCall;
            double HestonC = S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);

            // The output structure
            OutputMD output = new OutputMD();

            // Output the option price
            if(PutCall == "C")
                output.Price = HestonC;
            else
                output.Price = HestonP;

            // Output the integration domain;
            output.lower = A[0];
            output.upper = A[nj];

            // Output the number of integration points
            output.Npoints = (nj+1)*nX;
            return output;
        }
    }
}
