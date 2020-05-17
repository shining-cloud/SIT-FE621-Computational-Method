using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Heston
{
    class BisectionImpliedVol
    {
        // Bisection Algorithm for Black Scholes Implied Volatility =================================================================================
        public double BisectionMSIV(double S,double K,double r,double q,double T,double a,double b,double MktPrice,double Tol,int MaxIter,double B,double dt)
        {
            MSPrices MS = new MSPrices();
            double lowPrice = MS.MSPriceBS(S,K,T,a,r,q,MaxIter,Tol,B,dt);
            double highPrice = MS.MSPriceBS(S,K,T,b,r,q,MaxIter,Tol,B,dt);
            double lowCdif  = MktPrice - lowPrice;
            double highCdif = MktPrice - highPrice;
            double y = 0.0;
            double midP;
            if(lowCdif*highCdif > 0.0)
                y = -1.0;
            else
            {
                for(int x=0;x<=MaxIter;x++)
                {
                    midP = (a+b)/2.0;
                    double midCdif  = MktPrice - MS.MSPriceBS(S,K,T,midP,r,q,MaxIter,Tol,B,dt);
                    if(Math.Abs(midCdif) < Tol)
                    {
                        break;
                    }
                    else
                    {
                        if(midCdif > 0.0)
                            a = midP;
                        else
                            b = midP;

                    }
                    y = midP;
                }
            }
            return y;
        }
    }
}
