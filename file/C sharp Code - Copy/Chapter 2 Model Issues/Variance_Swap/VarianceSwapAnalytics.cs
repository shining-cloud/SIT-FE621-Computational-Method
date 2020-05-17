using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Variance_Swap
{
    class VarSwap
    {
        public double VarianceSwap(double[] KCI,double[] CallVI,double[] KPI,double[] PutVI,double S,double T,double rf,double q)
        {
            // Demeterfi et al (1999) fair strike of a variance swap
            // INPUTS
            //   KCI = Grid of call strikes
            //   CallVI = Interpolated call implied vol along the grid KCI
            //   KPI  = Grid of put strikes
            //   PutVI  = Interpolated put implied vol along the grid KPI
            //   S = Spot price
            //   T = Maturity
            //   rf = interest rate
            //   q = dividend yield

            // Do the required calculations on calls.  
            // Rename CallK and CallV for convenience.
            int n = CallVI.Length;

            // Take ATM as the boundary point
            double Sb = S;

            // Start the replication algorithm
            double[] Temp = new double[n-1];
            double[] CallWeight  = new double[n-1];
            double[] CallValue   = new double[n-1];
            double[] CallContrib = new double[n-1];
            BlackScholesPrice BS = new BlackScholesPrice();
            for(int k=0;k<=n-2;k++)
            {
                Temp[k] = (f(KCI[k+1],Sb,T) - f(KCI[k],Sb,T)) / (KCI[k+1] - KCI[k]);
                if(k==0)
                    CallWeight[0] = Temp[0];
                CallValue[k] = BS.BlackScholes(S,KCI[k],T,rf,q,CallVI[k],"C");
                if(k>0)
                    CallWeight[k] = Temp[k] - Temp[k-1];
                CallContrib[k] = CallValue[k]*CallWeight[k];
            }
            double Pi1 = CallContrib.Sum();

            // Do the calculations on puts. Flip the Vectors for Convenience
            n = PutVI.Length;
            Array.Reverse(KPI);
            Array.Reverse(PutVI);
            double[] Temp2 = new double[n-1];
            double[] PutWeight  = new double[n-1];
            double[] PutValue   = new double[n-1];
            double[] PutContrib = new double[n-1];
            for(int k=0;k<=n-2;k++)
            {
                Temp2[k] = (f(KPI[k+1],Sb,T) - f(KPI[k],Sb,T)) / (KPI[k] - KPI[k+1]);
                if(k==0)
                    PutWeight[0] = Temp2[0];
                PutValue[k] = BS.BlackScholes(S,KPI[k],T,rf,q,PutVI[k],"P");
                if(k>0)
                    PutWeight[k] = Temp2[k] - Temp2[k-1];
                PutContrib[k] = PutValue[k] * PutWeight[k];
            }
            double Pi2 = PutContrib.Sum();

            // Total cost of the portfolio
            double Pi_CP = Pi1 + Pi2;

            // Results of the replication
            // Estimate of fair variance
            double Kvar = 2.0/T*(rf*T - (S/Sb*Math.Exp(rf*T) - 1.0) - Math.Log(Sb/S)) + Math.Exp(rf*T)*Pi_CP;

            return Kvar;
        }

        // Auxiliary function
        public double f(double S,double Sb,double T)
        {
            return 2.0/T*((S - Sb) / Sb - Math.Log(S/Sb));
        }
    }
}


