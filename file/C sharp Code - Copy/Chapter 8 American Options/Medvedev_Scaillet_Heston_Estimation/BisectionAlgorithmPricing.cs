using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Heston
{ 
    class BisectionPricing
    {
        // Derivative of the MS approximation
        public double MSPutHestonDiff(double y,double theta,double Strike,HParam param,double r,double q,double T,int NumTerms,double dy)
        {
            MSExpansionHeston MS = new MSExpansionHeston();
            return (MS.MSPutHeston(y+dy,theta,Strike,param,r,q,T,NumTerms) - MS.MSPutHeston(y-dy,theta,Strike,param,r,q,T,NumTerms))/2.0/dy;
        }
        // Bisection Algorithm for MS Heston Prices
        public double BisectionMS(double a,double b,double theta,double Strike,HParam param,double r,double q,double T,int NumTerms,double Tol,int MaxIter,double dy)
        {
            double lowCdif  = MSPutHestonDiff(a,theta,Strike,param,r,q,T,NumTerms,dy);
            double highCdif = MSPutHestonDiff(b,theta,Strike,param,r,q,T,NumTerms,dy);
            double y = 0.0;
            double midP;
            if(lowCdif*highCdif > 0.0)
                y = -1.0;
            else
            {
                for(int x=0;x<=MaxIter;x++)
                {
                    midP = (a+b)/2.0;
                    double midCdif  = MSPutHestonDiff(midP,theta,Strike,param,r,q,T,NumTerms,dy);
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
        // Derivative of the MS approximation
        public double MSPutBSdiff(double y,MSsetBS mssettings,double dy)
        {
            MSExpansionBS MS = new MSExpansionBS();
            return (MS.MSPutBS(y+dy,mssettings) - MS.MSPutBS(y-dy,mssettings))/2.0/dy;
        }

        // Bisection Algorithm for Black Scholes pricing 
        public double BisectionBS(MSsetBS mssettings,double a,double b,double Tol,int MaxIter,double dy)
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
        // Standard Normal CDF ===================================================================================
        public double NormCDF(double x)
        {
            double x1 = 7.0*Math.Exp(-0.5*x*x);
            double x2 = 16.0*Math.Exp(-x*x*(2.0 - Math.Sqrt(2.0)));
            double x3 = (7.0 + 0.25*Math.PI*x*x)*Math.Exp(-x*x);
            double Q = 0.5*Math.Sqrt(1.0 - (x1 + x2 + x3)/30.0);
            if(x > 0)
                return 0.5 + Q;
            else
                return 0.5 - Q;
        }

        // Standard Normal PDF
        public double NormPDF(double x)
        {
            return Math.Exp(-0.5*x*x)/Math.Sqrt(2.0*Math.PI);
        }
    }
}
