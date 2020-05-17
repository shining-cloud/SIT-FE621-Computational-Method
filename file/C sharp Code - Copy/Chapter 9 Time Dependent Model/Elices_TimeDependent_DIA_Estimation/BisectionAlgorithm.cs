using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Elices_TimeDependent_DIA_Estimation
{
    class BisectionAlgo
    {
        // Bisection Algorithm for Black Scholes Implied Volatility
        public double BisecBSIV(string PutCall,double S,double K,double rf,double q,double T,double a,double b,double MktPrice,double Tol,int MaxIter)
        {
            BlackScholesPrice BS = new BlackScholesPrice();
            double lowCdif  = MktPrice - BS.BlackScholes(S,K,T,rf,q,a,PutCall);
	        double highCdif = MktPrice - BS.BlackScholes(S,K,T,rf,q,b,PutCall);
            double BSIV = 0.0;
            double midP;
            if (lowCdif*highCdif > 0.0)
	            BSIV = -1.0;
            else
            {
                for (int x=0; x<=MaxIter; x++)
                {
        		    midP = (a+b)/2.0;
                    double midCdif  = MktPrice - BS.BlackScholes(S,K,T,rf,q,midP,PutCall);
            		if (Math.Abs(midCdif) < Tol)
                    {
                		break;
                    }
                    else
                    {
            			if (midCdif > 0.0)
            				a = midP;
	            		else
            				b = midP;
                    }
                    BSIV = midP;
                }
            }
            return BSIV;
        }
        // Random number in (a,b) ==========================================================
        public readonly Random U = new Random();
        public readonly object sync = new object();
        public double RandomNum(double a,double b)
        {
            int divisor = 1000000000;
            lock(sync) { return a + (b-a)*U.Next(0,divisor)/divisor; }
        }
    }
}
