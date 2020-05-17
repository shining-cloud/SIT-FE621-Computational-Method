using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Black_Scholes
{
    class BisectionAlgo
    {
        // Derivative of the MS approximation
        public double MSPutBSdiff(double y,MSset mssettings,double dy)
        {
            MSPrice MS = new MSPrice();
            return (MS.MSPutBS(y+dy,mssettings) - MS.MSPutBS(y-dy,mssettings))/2.0/dy;
        }

        // Bisection Algorithm for Black Scholes pricing 
        public double Bisection(MSset mssettings,double a,double b,double Tol,int MaxIter,double dy)
        {
            double lowCdif  = MSPutBSdiff(a,mssettings,dy);
            double highCdif = MSPutBSdiff(b,mssettings,dy);
            double y = 0.0;
            double midP;
            if(lowCdif*highCdif > 0.0)
                y = -999.0;
            else
            {
                for(int x=0;x<=MaxIter;x++)
                {
                    midP = (a+b)/2.0;
                    double midCdif  = MSPutBSdiff(midP,mssettings,dy);
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
