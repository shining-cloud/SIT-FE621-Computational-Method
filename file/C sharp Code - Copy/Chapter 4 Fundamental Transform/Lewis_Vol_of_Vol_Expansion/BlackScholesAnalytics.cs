using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Vol_of_Vol_Expansion
{
    class BlackScholes
    {
        // Black Scholes Call, using variance as input
        public double BSC(double S,double K,double rf,double q,double v,double T)
        {
            double d1 = (Math.Log(S/K) + T*(rf - q + 0.5*v)) / Math.Sqrt(v*T);
            double d2 = d1 - Math.Sqrt(v*T);
            return S*Math.Exp(-q*T)*NormCDF(d1) - Math.Exp(-rf*T)*K*NormCDF(d2);
        }
        // Black Scholes "vega" -- derivatives w.r.t. variance v
        public double BSV(double S,double K,double rf,double q,double v,double T)
        {
            double pi = Math.PI;
            return Math.Sqrt(T/8.0/pi/v)*S*Math.Exp(-q*T)*Math.Exp(-0.5*Math.Pow(((Math.Log(S/K) + (rf-q+v/2.0)*T)/Math.Sqrt(v*T)),2));
        }
        // Standard Normal CDF
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
    }
}
