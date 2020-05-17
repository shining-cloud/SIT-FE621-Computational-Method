using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Heston
{
    class GoldenSearch
    {
        // Golden Section Search method for Black Scholes
        public double GoldenSearchBS(double a,double b,double tol,int MaxIter,MSsetBS mssettings)
        {
            MSExpansionBS MS = new MSExpansionBS();
            double x1,x2,x3,x4,f1,f2,f3,f4;
            double GR = (Math.Sqrt(5.0) - 1.0)/2.0;
            int k = 0;
            while(Math.Abs(b-a) > tol)
            {
                k = k+1;
                x1 = a;
                x2 = a + (1.0-GR)*(b-a);
                x3 = a + GR*(b-a);
                x4 = b;
                f1 = -MS.MSPutBS(x1,mssettings);
                f2 = -MS.MSPutBS(x2,mssettings);
                f3 = -MS.MSPutBS(x3,mssettings);
                f4 = -MS.MSPutBS(x4,mssettings);
                if((f1 > f2) & (f2 < f3))
                    b = x3;
                else if((f2 > f3) & (f3 < f4))
                    a = x2;
                if(k > MaxIter)
                    break;
            }
            return (a+b)/2.0;
        }
        // Golden Section search form Heston
        public double GoldenSearchMS(double a,double b,double tol, int MaxIter, double K,HParam param,double theta,double r,double q,double T,int NumTerms)
        {
            MSExpansionHeston MS = new MSExpansionHeston();
            double x1,x2,x3,x4,f1,f2,f3,f4;
            double GR = (Math.Sqrt(5.0) - 1.0)/2.0;
            int k = 0;
            while(Math.Abs(b-a) > tol)
            {
                k = k+1;
                x1 = a;
                x2 = a + (1.0-GR)*(b-a);
                x3 = a + GR*(b-a);
                x4 = b;
                f1 = -MS.MSPutHeston(x1,theta,K,param,r,q,T,NumTerms);
                f2 = -MS.MSPutHeston(x2,theta,K,param,r,q,T,NumTerms);
                f3 = -MS.MSPutHeston(x3,theta,K,param,r,q,T,NumTerms);
                f4 = -MS.MSPutHeston(x4,theta,K,param,r,q,T,NumTerms);
                if((f1 > f2) & (f2 < f3))
                    b = x3;
                else if((f2 > f3) & (f3 < f4))
                    a = x2;
                if(k > MaxIter)
                    break;
            }
            return (a+b)/2.0;
        }
    }
}



