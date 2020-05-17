using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Gauthier_Starting_Values
{
    class BlackScholesPrice
    {
        // Standard Normal CDF ==============================================================
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
        // Black Scholes Price of Call or put ====================================================================
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
        // Black Scholes derivatives
        public double[] BlackScholesDerivatives(double kappa,double theta,double v0,double S,double K,double T,double rf,double q)
        {
            // log Spot price
            double x = Math.Log(S);

            // Integrated variance
            double wT = (v0-theta)*(1-Math.Exp(-kappa*T))/kappa + theta*T;
            double y = wT;

            // Black Scholes Put Price
            double g = Math.Pow(y,-0.5)*(-x + Math.Log(K) - (rf-q)*T) - 0.5*Math.Sqrt(y);
            double f = Math.Pow(y,-0.5)*(-x + Math.Log(K) - (rf-q)*T) + 0.5*Math.Sqrt(y);
            double BSPut = K*Math.Exp(-rf*T)*NormCDF(f) - S*Math.Exp(-q*T)*NormCDF(g);

            // Normal pdf, phi(f) and phi(g)

            double pi = Math.PI;
            double phif = Math.Exp(-f*f/2.0)/Math.Sqrt(2.0*pi);
            double phig = Math.Exp(-g*g/2.0)/Math.Sqrt(2.0*pi);

            // Derivatives of f and g
            double fx = -Math.Pow(y,-0.5);
            double fy = -1.0/2.0/y*g;
            double gx = fx;
            double gy = -1.0/2.0/y*f;

            // The cdf PHI(f) and PHI(g)
            double PHIf = NormCDF(f);
            double PHIg = NormCDF(g);

            // Derivatives of the pdf phi(f)
            double phifx = Math.Pow(y,-0.5)*f*phif;
            double phify = 0.5/y*f*g*phif;

            // Derivatives of the cdf PHI(f)
            double PHIfxy   = 0.5*Math.Pow(y,-1.5)*phif*(1-f*g);
            double PHIfx2y  = 0.5*Math.Pow(y,-2.0)*phif*(2*f+g-f*f*g);
            double PHIfy2   = 0.5*Math.Pow(y,-2.0)*phif*(g+f/2.0-f*g*g/2);
            double PHIfx2y2 = 0.5*((Math.Pow(y,-2.0)*phify-2.0*Math.Pow(y,-3.0)*phif)*(2.0*f+g-f*f*g) + Math.Pow(y,-2.0)*phif*(2.0*fy+gy-2.0*f*fy*g-f*f*gy));

            // Derivatives of the pdf phi(g)
            double phigx = Math.Pow(y,-0.5)*g*phig;
            double phigy = 0.5/y*f*g*phig;

            // Derivatives of cdf PHI(g)
            double PHIgx = -phig*Math.Pow(y,-0.5);
            double PHIgy   = -0.5*f*phig/y;
            double PHIgxy  =  0.5*Math.Pow(y,-1.5)*phig*(1.0-f*g);
            double PHIgx2y =  0.5*Math.Pow(y,-2.0)*phig*(2.0*g+f-g*g*f);
            double PHIgy2  =  0.5*Math.Pow(y,-2.0)*phig*(f+g/2.0-g*f*f/2);
            double PHIgxy2 =  0.5*Math.Pow(y,-2.0)*(phigx*(f+g/2.0-f*f*g/2.0) + phig*(fx+gx/2.0-f*fx*g/2.0-f*f*gx/2.0));
            double PHIgx2y2 = 0.5*((Math.Pow(y,-2.0)*phigy-2*Math.Pow(y,-3.0)*phig)*(2.0*g+f-g*g*f) + Math.Pow(y,-2.0)*phig*(2.0*gy+fy-2.0*g*gy*f-g*g*fy));

            // Derivatives of Black-Scholes Put
            double dPdxdy   = K*Math.Exp(-rf*T)*PHIfxy   - Math.Exp(-q*T)*S*(PHIgy + PHIgxy);
            double dPdx2dy  = K*Math.Exp(-rf*T)*PHIfx2y  - Math.Exp(-q*T)*S*(PHIgy + 2*PHIgxy + PHIgx2y);
            double dPdy2    = K*Math.Exp(-rf*T)*PHIfy2   - Math.Exp(-q*T)*S*PHIgy2;
            double dPdx2dy2 = K*Math.Exp(-rf*T)*PHIfx2y2 - Math.Exp(-q*T)*S*(PHIgy + 2*PHIgxy + PHIgx2y + PHIgy2 + 2*PHIgxy2 + PHIgx2y2);

            // Return the quantities
            double[] output = new double[4];
            output[0] = dPdxdy;
            output[1] = dPdx2dy;
            output[2] = dPdy2;
            output[3] = dPdx2dy2;

            return output;
        }
    }
}