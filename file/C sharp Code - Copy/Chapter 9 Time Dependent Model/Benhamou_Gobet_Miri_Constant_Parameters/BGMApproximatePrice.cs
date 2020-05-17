using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Benhamou_Gobet_Miri_Constant_Parameters
{
    class BGMPrice
    {
        public double BGMApproxPrice(HParam param,OpSet settings,double K,double T)
        {
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = 0.0;

            double S = settings.S;
            double rf = settings.r;
            double q = settings.q;
            string PutCall = settings.PutCall;

            // Log Spot price
            double x = Math.Log(settings.S);

            // Integrated variance
            double wT = (v0-theta)*(1-Math.Exp(-kappa*T))/kappa + theta*T;
            double y = wT;

            // Black Scholes Put Price
            BlackScholesPrice BS = new BlackScholesPrice();
            double g = Math.Pow(y,-0.5) * (-x + Math.Log(K) - (rf-q)*T) - 0.5*Math.Sqrt(y);
            double f = Math.Pow(y,-0.5) * (-x + Math.Log(K) - (rf-q)*T) + 0.5*Math.Sqrt(y);
            double BSPut = K*Math.Exp(-rf*T)*BS.NormCDF(f) - S*Math.Exp(-q*T)*BS.NormCDF(g);

            // Shortcut notation
            double k  = kappa;
            double kT = kappa*T;
            double ekT  = Math.Exp(k*T);
            double ekTm = Math.Exp(-k*T);

            // Coefficients for the expansion
            double a1T = (rho*sigma*ekTm/k/k) * (v0*(-kT+ekT-1.0) + theta*(kT+ekT*(kT-2.0)+2.0));
            double a2T = (rho*rho*sigma*sigma*ekTm/2.0/(k*k*k)) * (v0*(-kT*(kT+2.0)+2.0*ekT-2.0) + theta*(2.0*ekT*(kT-3.0)+kT*(kT+4.0)+6.0));
            double b0T = (sigma*sigma*Math.Exp(-2.0*kT)/4.0/(k*k*k)) * (v0*(-4.0*ekT*kT+2.0*Math.Exp(2.0*kT)-2.0) + theta*(4.0*ekT*(kT+1.0)+Math.Exp(2.0*kT)*(2.0*kT-5.0)+1.0));
            double b2T = a1T*a1T/2.0;

            // Normal pdf, phi(f) and phi(g)
            double pi = Math.PI;
            double phif = Math.Exp(-f*f/2.0)/Math.Sqrt(2.0*pi);
            double phig = Math.Exp(-g*g/2.0)/Math.Sqrt(2.0*pi);

            // Derivatives of f and g
            double fx = -Math.Pow(y,-0.5);
            double fy = -0.5/y*g;
            double gx = fx;
            double gy = -0.5/y*f;

            // The cdf PHI(f) and PHI(g)
            double PHIf = BS.NormCDF(f);
            double PHIg = BS.NormCDF(g);

            // Derivatives of the pdf phi(f)
            double phifx = Math.Pow(y,-0.5)*f*phif;
            double phify = 0.5/y*f*g*phif;

            // Derivatives of the cdf PHI(f)
            double PHIfxy   = 0.5*Math.Pow(y,-1.5)*phif*(1.0-f*g);
            double PHIfx2y  = 0.5*Math.Pow(y,-2.0)*phif*(2.0*f+g-f*f*g);
            double PHIfy2   = 0.5*Math.Pow(y,-2.0)*phif*(g+f/2.0-f*g*g/2.0);
            double PHIfx2y2 = 0.5*((Math.Pow(y,-2.0)*phify-2.0*Math.Pow(y,-3.0)*phif)*(2.0*f+g-f*f*g) +
	               Math.Pow(y,-2.0)*phif*(2.0*fy+gy-2.0*f*fy*g-f*f*gy));

            // Derivatives of the pdf phi(g)
            double phigx = Math.Pow(y,-0.5)*g*phig;
            double phigy = 0.5/y*f*g*phig;

            // Derivatives of cdf PHI(g)
            double PHIgx = -phig*Math.Pow(y,-0.5);
            double PHIgy   = -0.5*f*phig/y;
            double PHIgxy  =  0.5*Math.Pow(y,-1.5)*phig*(1.0-f*g);
            double PHIgx2y =  0.5*Math.Pow(y,-2.0)*phig*(2.0*g+f-g*g*f);
            double PHIgy2  =  0.5*Math.Pow(y,-2.0)*phig*(f+g/2.0-g*f*f/2.0);
            double PHIgxy2 =  0.5*Math.Pow(y,-2.0)*(phigx*(f+g/2.0-f*f*g/2.0) + phig*(fx+gx/2.0-f*fx*g/2.0-f*f*gx/2.0));
            double PHIgx2y2 = 0.5*((Math.Pow(y,-2.0)*phigy-2.0*Math.Pow(y,-3.0)*phig)*(2.0*g+f-g*g*f) +
                   Math.Pow(y,-2.0)*phig*(2.0*gy+fy-2.0*g*gy*f-g*g*fy));

            // Derivatives of Black-Scholes Put
            double dPdxdy   = K*Math.Exp(-rf*T)*PHIfxy   - Math.Exp(-q*T)*S*(PHIgy + PHIgxy);
            double dPdx2dy  = K*Math.Exp(-rf*T)*PHIfx2y  - Math.Exp(-q*T)*S*(PHIgy + 2.0*PHIgxy + PHIgx2y);
            double dPdy2    = K*Math.Exp(-rf*T)*PHIfy2   - Math.Exp(-q*T)*S*PHIgy2;
            double dPdx2dy2 = K*Math.Exp(-rf*T)*PHIfx2y2 - Math.Exp(-q*T)*S*(PHIgy + 2.0*PHIgxy + PHIgx2y + PHIgy2 + 2.0*PHIgxy2 + PHIgx2y2);

            // Benhamou, Gobet, Miri expansion
            double Put = BSPut + a1T*dPdxdy + a2T*dPdx2dy + b0T*dPdy2 + b2T*dPdx2dy2;

            // Return the put or the call by put-call parity
            if(PutCall == "P")
                return Put;
            else
                return Put - K*Math.Exp(-rf*T) + S*Math.Exp(-q*T);
        }
    }
}

	
