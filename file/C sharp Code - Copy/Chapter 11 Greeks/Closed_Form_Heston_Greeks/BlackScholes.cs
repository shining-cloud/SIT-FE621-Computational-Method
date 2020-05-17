using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Closed_Form_Heston_Greeks
{
    class BlackScholesPrice
    {
        // Standard Normal CDF
        public double NormCDF(double x)
        {
            double b0 =  0.2316419;
            double b1 =  0.319381530;
            double b2 = -0.356563782;
            double b3 =  1.781477937;
            double b4 = -1.821255978;
            double b5 =  1.330274429;
            double pi = Math.PI;
            double phi = Math.Exp(-x*x/2.0)/Math.Sqrt(2.0*pi);
            double t,c;
            double CDF = 0.5;
            if(x > 0.0)
            {
                t = 1.0/(1.0 + b0*x);
                CDF = 1.0 - phi*(b1*t + b2*Math.Pow(t,2) + b3*Math.Pow(t,3) + b4*Math.Pow(t,4) + b5*Math.Pow(t,5));
            }
            else if(x < 0.0)
            {
                x = -x;
                t = 1.0/(1.0 + b0*x);
                c = 1.0 - phi*(b1*t + b2*Math.Pow(t,2) + b3*Math.Pow(t,3) + b4*Math.Pow(t,4) + b5*Math.Pow(t,5));
                CDF = 1.0 - c;
            }
            return CDF;
        }
        // Standard Normal PDF
        public double NormPDF(double x)
        {
            return Math.Exp(-x*x/2.0) / Math.Sqrt(2*Math.PI);
        }
        // Black Scholes Price of Call or put
        public double BlackScholes(double S,double K,double T,double rf,double q,double v,string PutCall)
        {
            double d1 = (Math.Log(S/K) + (rf-q+v*v/2.0)*T) / v / Math.Sqrt(T);
            double d2 = d1 - v*Math.Sqrt(T);
            double BSCall = S*Math.Exp(-q*T)*NormCDF(d1) - K*Math.Exp(-rf*T)*NormCDF(d2);
            double Price = 0.0;
            if(PutCall == "C")
            {
                Price = BSCall;
            }
            else if(PutCall == "P")
            {
                Price = BSCall - S*Math.Exp(-q*T) + K*Math.Exp(-rf*T);
            }
            return Price;
        }
        // Black Scholes Greeks
        public double BSGreeks(string PutCall,double S,double K,double r,double q,double T,double sigma,string Greek)
        {
            double d1 = (Math.Log(S/K) + (r-q+sigma*sigma/2.0)*T)/sigma/Math.Sqrt(T);
            double d2 = d1 - sigma*Math.Sqrt(T);
            double Nd1 = NormCDF(d1);
            double Nd2 = NormCDF(d2);
            double nd1 = NormPDF(d1);
            double y = 0.0;

            if (Greek == "Delta")
            {
                if(PutCall == "C")
                    y = Math.Exp(-q*T)*Nd1;
                else
                    y = Math.Exp(-q*T)*(Nd1 - 1.0);
            }
            else if(Greek == "Gamma")
                y = nd1*Math.Exp(-q*T)/S/sigma/Math.Sqrt(T);
            else if(Greek == "Vega")
                y = S*Math.Exp(-q*T)*Math.Sqrt(T)*nd1;
            else if(Greek == "Rho")
            {
                if(PutCall == "C")
                    y = K*T*Math.Exp(-r*T)*Nd2;
                else
                    y = K*T*Math.Exp(-r*T)*(Nd2-1);
            }
            else if(Greek == "Theta")
            {
                if(PutCall == "C")
                    y = -Math.Exp(-q*T)*S*nd1*sigma/2.0/Math.Sqrt(T) - r*K*Math.Exp(-r*T)*Nd2 + q*S*Math.Exp(-q*T)*Nd1;
                else
                    y = -Math.Exp(-q*T)*S*nd1*sigma/2.0/Math.Sqrt(T) + r*K*Math.Exp(-r*T)*(1-Nd2) - q*S*Math.Exp(-q*T)*(1-Nd1);
            }
            else if(Greek == "Vanna")
                y = -Math.Exp(-q*T)*nd1*d2/sigma;
            else if(Greek == "Volga")
                y = S*Math.Exp(-q*T)*nd1*Math.Sqrt(T)*d1*d2/sigma;
            return y;
        }
    }
}

