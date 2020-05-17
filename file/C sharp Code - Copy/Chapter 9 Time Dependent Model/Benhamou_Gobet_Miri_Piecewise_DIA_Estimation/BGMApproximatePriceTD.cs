using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Benhamou_Gobet_Miri_Piecewise_DIA_Estimation
{
    class BGMPrice
    {
        // Parameters for the shortest maturity appear first in the vectors theta, sigma, and rho.
        // Parameters for the longest maturity appear last.
        public double BGMApproxPriceTD(double kappa,double v0,double[] theta,double[] sigma,double[] rho,OPSet settings,double K,double[] T)
        {
            int NT = T.Length;
            double[] a1T = new double[NT];
            double[] a2T = new double[NT];
            double[] b0T = new double[NT];
            double[] b2T = new double[NT];

            // First set of coefficients
            a1T[0] = -rho[0]*sigma[0]*(2.0*theta[0]*Math.Exp(kappa*T[0])+v0*T[0]*kappa+v0-theta[0]*T[0]*kappa*Math.Exp(kappa*T[0])-theta[0]*kappa*T[0]-2.0*theta[0]-v0*Math.Exp(kappa*T[0]))/Math.Exp(kappa*T[0])/Math.Pow(kappa,2.0);
            a2T[0] = -0.5*rho[0]*rho[0]*sigma[0]*sigma[0]*(Math.Pow(kappa,2.0)*v0*Math.Pow(T[0],2.0)+6.0*theta[0]*Math.Exp(kappa*T[0])+2.0*v0*T[0]*kappa+2.0*v0-Math.Pow(kappa,2.0)*theta[0]*Math.Pow(T[0],2.0)-4.0*theta[0]*kappa*T[0]-2.0*theta[0]*T[0]*kappa*Math.Exp(kappa*T[0])-6.0*theta[0]-2.0*v0*Math.Exp(kappa*T[0]))/Math.Exp(kappa*T[0])/Math.Pow(kappa,3.0);
            b0T[0] = -0.25*sigma[0]*sigma[0]*(5.0*theta[0]*Math.Pow(Math.Exp(kappa*T[0]), 2.0)+4.0*Math.Exp(kappa*T[0])*v0*T[0]*kappa-4.0*theta[0]*T[0]*kappa*Math.Exp(kappa*T[0])-2.0*theta[0]*T[0]*kappa*Math.Pow(Math.Exp(kappa*T[0]), 2.0)+2.0*v0-theta[0]-2.0*v0*Math.Pow(Math.Exp(kappa*T[0]), 2.0)-4.0*theta[0]*Math.Exp(kappa*T[0]))/Math.Pow(Math.Exp(kappa*T[0]), 2.0)/Math.Pow(kappa,3.0);
            b2T[0] =  a1T[0]*a1T[0]/2.0;

            double A1=0.0,A2=0.0,B0=0.0;
            int j;
            // Coefficients under multiple maturities
            if(NT>=2)
            {
                for(int t=0;t<=NT-2;t++)
                {
                    // Coefficients for a1T[t+1]
                    j = 0;
                        A1 += rho[j]*sigma[j]*(-Math.Exp(-kappa*T[t])+Math.Exp(-kappa*T[t+1]))*(theta[j]-v0*T[j]*kappa+theta[j]*T[j]*kappa-theta[j]*Math.Exp(kappa*T[j]))/Math.Pow(kappa,2.0);
                    for(j=1;j<=t;j++)
                        A1 += -rho[j]*sigma[j]*(-Math.Exp(-kappa*T[t])+Math.Exp(-kappa*T[t+1]))*(-v0*T[j-1]*kappa+theta[j]*T[j-1]*kappa-theta[j]*Math.Exp(kappa*T[j-1])+v0*T[j]*kappa-theta[j]*T[j]*kappa+theta[j]*Math.Exp(kappa*T[j]))/Math.Pow(kappa,2.0);
                    j = t+1;
                        A1 += -rho[j]*sigma[j]*(-v0*Math.Exp(kappa*T[t+1])-v0*T[t]*kappa*Math.Exp(kappa*T[t])+theta[j]*Math.Exp(kappa*T[t+1])+theta[j]*T[t]*kappa*Math.Exp(kappa*T[t])+theta[j]*T[t]*kappa*Math.Exp(kappa*(T[t]+T[t+1]))-theta[j]*Math.Exp(2.0*kappa*T[t])+Math.Exp(kappa*T[t])*v0-Math.Exp(kappa*(T[t]+T[t+1]))*theta[j]*T[t+1]*kappa+theta[j]*Math.Exp(kappa*(T[t]+T[t+1]))+Math.Exp(kappa*T[t])*v0*kappa*T[t+1]-theta[j]*Math.Exp(kappa*T[t])-Math.Exp(kappa*T[t])*theta[j]*kappa*T[t+1])/Math.Pow(kappa,2.0)*Math.Exp(-kappa*(T[t]+T[t+1]));
                    a1T[t+1] = a1T[t] + A1;

                    // Coefficients for a2T[t+1]
                    j = 0;
                        A2 += -Math.Pow(rho[j],2.0)*Math.Pow(sigma[j],2.0)*(Math.Exp(-kappa*T[t])+Math.Exp(-kappa*T[t+1])*T[t]*kappa-Math.Exp(-kappa*T[t+1])-Math.Exp(-kappa*T[t+1])*kappa*T[t+1])*(theta[j]-v0*T[j]*kappa+theta[j]*T[j]*kappa-theta[j]*Math.Exp(kappa*T[j]))/Math.Pow(kappa,3.0)
                              + 1.0/2.0*Math.Pow(rho[j],2.0)*Math.Pow(sigma[j],2.0)*(-Math.Exp(-kappa*T[t])+Math.Exp(-kappa*T[t+1]))*(2.0*kappa*theta[j]*T[t]+2.0*theta[j]-2.0*v0*T[t]*T[j]*Math.Pow(kappa,2.0)+v0*Math.Pow(T[j],2.0)*Math.Pow(kappa,2.0)+2.0*theta[j]*T[t]*T[j]*Math.Pow(kappa,2.0)-theta[j]*Math.Pow(T[j],2.0)*Math.Pow(kappa,2.0)-2.0*theta[j]*T[t]*Math.Exp(kappa*T[j])*kappa+2.0*theta[j]*Math.Exp(kappa*T[j])*kappa*T[j]-2.0*theta[j]*Math.Exp(kappa*T[j]))/Math.Pow(kappa,3.0);
                    for(j=1;j<=t;j++)
                        A2 += Math.Pow(rho[j],2.0)*Math.Pow(sigma[j],2.0)*(Math.Exp(-kappa*T[t])+Math.Exp(-kappa*T[t+1])*T[t]*kappa-Math.Exp(-kappa*T[t+1])-Math.Exp(-kappa*T[t+1])*kappa*T[t+1])*(-v0*T[j-1]*kappa+theta[j]*T[j-1]*kappa-theta[j]*Math.Exp(kappa*T[j-1])+v0*T[j]*kappa-theta[j]*T[j]*kappa+theta[j]*Math.Exp(kappa*T[j]))/Math.Pow(kappa,3.0)
                            - 1.0/2.0*Math.Pow(rho[j],2.0)*Math.Pow(sigma[j],2.0)*(-Math.Exp(-kappa*T[t])+Math.Exp(-kappa*T[t+1]))*(-2.0*v0*T[t]*T[j-1]*Math.Pow(kappa,2.0)+v0*Math.Pow(T[j-1],2.0)*Math.Pow(kappa,2.0)+2.0*theta[j]*T[t]*T[j-1]*Math.Pow(kappa,2.0)-theta[j]*Math.Pow(T[j-1],2.0)*Math.Pow(kappa,2.0)-2.0*theta[j]*T[t]*Math.Exp(kappa*T[j-1])*kappa+2.0*theta[j]*Math.Exp(kappa*T[j-1])*kappa*T[j-1]-2.0*theta[j]*Math.Exp(kappa*T[j-1])+2.0*v0*T[t]*T[j]*Math.Pow(kappa,2.0)-v0*Math.Pow(T[j],2.0)*Math.Pow(kappa,2.0)-2.0*theta[j]*T[t]*T[j]*Math.Pow(kappa,2.0)+theta[j]*Math.Pow(T[j],2.0)*Math.Pow(kappa,2.0)+2.0*theta[j]*T[t]*Math.Exp(kappa*T[j])*kappa-2.0*theta[j]*Math.Exp(kappa*T[j])*kappa*T[j]+2.0*theta[j]*Math.Exp(kappa*T[j]))/Math.Pow(kappa,3.0);
                    j = t+1;
                        A2 += 1.0/2.0*Math.Pow(rho[j],2.0)*Math.Pow(sigma[j],2.0)*(2.0*v0*Math.Exp(kappa*T[t+1])-v0*Math.Pow(T[t],2.0)*Math.Pow(kappa,2.0)*Math.Exp(kappa*T[t])+2.0*v0*T[t]*kappa*Math.Exp(kappa*T[t])+2.0*v0*Math.Pow(kappa,2.0)*T[t+1]*T[t]*Math.Exp(kappa*T[t])-2.0*theta[j]*Math.Exp(kappa*T[t+1])+theta[j]*Math.Pow(T[t],2.0)*Math.Pow(kappa,2.0)*Math.Exp(kappa*T[t])-2.0*theta[j]*T[t]*kappa*Math.Exp(kappa*T[t])-2.0*theta[j]*Math.Pow(kappa,2.0)*T[t+1]*T[t]*Math.Exp(kappa*T[t])-2.0*theta[j]*T[t]*kappa*Math.Exp(kappa*(T[t]+T[t+1]))-2.0*theta[j]*Math.Exp(2.0*kappa*T[t])*kappa*T[t]+4.0*theta[j]*Math.Exp(2.0*kappa*T[t])+2.0*theta[j]*T[t+1]*Math.Exp(2.0*kappa*T[t])*kappa-2.0*Math.Exp(kappa*T[t])*v0-2.0*Math.Exp(kappa*T[t])*v0*kappa*T[t+1]+2.0*Math.Exp(kappa*(T[t]+T[t+1]))*theta[j]*T[t+1]*kappa-Math.Exp(kappa*T[t])*Math.Pow(kappa,2.0)*v0*Math.Pow(T[t+1],2.0)-4.0*theta[j]*Math.Exp(kappa*(T[t]+T[t+1]))+Math.Exp(kappa*T[t])*Math.Pow(kappa,2.0)*theta[j]*Math.Pow(T[t+1],2.0)+2.0*theta[j]*Math.Exp(kappa*T[t])+2.0*Math.Exp(kappa*T[t])*theta[j]*kappa*T[t+1])/Math.Pow(kappa,3.0)*Math.Exp(-kappa*(T[t]+T[t+1]));
                    a2T[t+1] = a2T[t] + A2;

                    // Coefficients for b0T[t+1]
                    j = 0;
                        B0 += -1.0/4.0*Math.Pow(sigma[j],2.0)*(-Math.Exp(-2.0*kappa*T[t])+2.0*Math.Exp(-kappa*(T[t]+T[t+1]))-Math.Exp(-2.0*kappa*T[t+1]))*(-2.0*v0+theta[j]+2.0*v0*Math.Exp(kappa*T[j])-2.0*theta[j]*Math.Exp(kappa*T[j])+theta[j]*Math.Exp(2.0*kappa*T[j]))/Math.Pow(kappa,3.0)
                                + 1.0/2.0*Math.Pow(sigma[j],2.0)*(-theta[j]*Math.Exp(kappa*T[t+1])+theta[j]*Math.Exp(kappa*T[t])+2.0*v0*Math.Exp(kappa*T[t+1])-2.0*Math.Exp(kappa*T[t])*v0-2.0*theta[j]*Math.Exp(kappa*(T[t]+T[t+1]))+2.0*theta[j]*Math.Exp(2.0*kappa*T[t])+2.0*Math.Exp(kappa*(T[t]+T[t+1]))*v0*T[j]*kappa-2.0*Math.Exp(kappa*(T[t+1]+T[j]))*v0-2.0*Math.Exp(kappa*(T[t]+T[t+1]))*theta[j]*T[j]*kappa+2.0*Math.Exp(kappa*(T[t+1]+T[j]))*theta[j]+2.0*Math.Exp(kappa*(T[t]+T[t+1]+T[j]))*theta[j]-Math.Exp(kappa*(T[t+1]+2.0*T[j]))*theta[j]-2.0*Math.Exp(2.0*kappa*T[t])*v0*T[j]*kappa+2.0*Math.Exp(kappa*(T[t]+T[j]))*v0+2.0*Math.Exp(2.0*kappa*T[t])*theta[j]*T[j]*kappa-2.0*Math.Exp(kappa*(T[t]+T[j]))*theta[j]-2.0*Math.Exp(kappa*(2.0*T[t]+T[j]))*theta[j]+Math.Exp(kappa*(T[t]+2.0*T[j]))*theta[j])*Math.Exp(-kappa*(2.0*T[t]+T[t+1]))/Math.Pow(kappa,3.0);
                    for(j=1;j<=t;j++)
                        B0 += -1.0/4.0*Math.Pow(sigma[j],2.0)*(-Math.Exp(-2.0*kappa*T[t])+2.0*Math.Exp(-kappa*(T[t]+T[t+1]))-Math.Exp(-2.0*kappa*T[t+1]))*(-2.0*v0*Math.Exp(kappa*T[j-1])+2.0*theta[j]*Math.Exp(kappa*T[j-1])-theta[j]*Math.Exp(2.0*kappa*T[j-1])+2.0*v0*Math.Exp(kappa*T[j])-2.0*theta[j]*Math.Exp(kappa*T[j])+theta[j]*Math.Exp(2.0*kappa*T[j]))/Math.Pow(kappa,3.0)
                              + 1.0/2.0*Math.Pow(sigma[j],2.0)*(Math.Exp(-kappa*T[t])-Math.Exp(-kappa*T[t+1]))*(-2.0*v0*T[j-1]*kappa*Math.Exp(kappa*T[t])+2.0*v0*Math.Exp(kappa*T[j-1])+2.0*theta[j]*T[j-1]*kappa*Math.Exp(kappa*T[t])-2.0*theta[j]*Math.Exp(kappa*T[j-1])-2.0*theta[j]*Math.Exp(kappa*(T[j-1]+T[t]))+theta[j]*Math.Exp(2.0*kappa*T[j-1])+2.0*v0*T[j]*kappa*Math.Exp(kappa*T[t])-2.0*v0*Math.Exp(kappa*T[j])-2.0*theta[j]*T[j]*kappa*Math.Exp(kappa*T[t])+2.0*theta[j]*Math.Exp(kappa*T[j])+2.0*Math.Exp(kappa*(T[t]+T[j]))*theta[j]-theta[j]*Math.Exp(2.0*kappa*T[j]))*Math.Exp(-kappa*T[t])/Math.Pow(kappa,3.0);
                    j = t+1;
                        B0 += -1.0/4.0*Math.Pow(sigma[j],2.0)*(-2.0*v0*Math.Exp(2.0*kappa*T[t+1])-4.0*v0*T[t]*kappa*Math.Exp(kappa*(T[t]+T[t+1]))+2.0*v0*Math.Exp(2.0*kappa*T[t])+2.0*theta[j]*Math.Exp(2.0*kappa*T[t+1])+4.0*theta[j]*T[t]*kappa*Math.Exp(kappa*(T[t]+T[t+1]))-2.0*theta[j]*Math.Exp(2.0*kappa*T[t])+2.0*theta[j]*T[t]*kappa*Math.Exp(kappa*(T[t]+2.0*T[t+1]))-4.0*theta[j]*Math.Exp(kappa*(2.0*T[t]+T[t+1]))+theta[j]*Math.Exp(3*kappa*T[t])+3.0*theta[j]*Math.Exp(kappa*(T[t]+2.0*T[t+1]))+4.0*Math.Exp(kappa*(T[t]+T[t+1]))*v0*kappa*T[t+1]-4.0*Math.Exp(kappa*(T[t]+T[t+1]))*theta[j]*T[t+1]*kappa-2.0*Math.Exp(kappa*(T[t]+2.0*T[t+1]))*theta[j]*T[t+1]*kappa)/Math.Pow(kappa,3.0)*Math.Exp(-kappa*(T[t]+2.0*T[t+1]));
                    b0T[t+1] = b0T[t] + B0;

                    // Coefficients for b2T[t+1]
                    b2T[t+1] = a1T[t+1]*a1T[t+1]/2.0;
                }
            }

            // Coefficients for the expansion are the last ones in the iterations
            double A1T = a1T[NT-1];
            double A2T = a2T[NT-1];
            double B0T = b0T[NT-1];
            double B2T = b2T[NT-1];

            double S = settings.S;
            double rf = settings.r;
            double q = settings.q;
            string PutCall = settings.PutCall;

            // Log Spot price
            double x = Math.Log(settings.S);

            // Integrated variance
            double wT = (v0-theta[NT-1])*(1-Math.Exp(-kappa*T[NT-1]))/kappa + theta[NT-1]*T[NT-1];
            double y = wT;

            // Black Scholes Put Price
            BlackScholesPrice BS = new BlackScholesPrice();
            double g = Math.Pow(y,-0.5) * (-x + Math.Log(K) - (rf-q)*T[NT-1]) - 0.5*Math.Sqrt(y);
            double f = Math.Pow(y,-0.5) * (-x + Math.Log(K) - (rf-q)*T[NT-1]) + 0.5*Math.Sqrt(y);
            double BSPut = K*Math.Exp(-rf*T[NT-1])*BS.NormCDF(f) - S*Math.Exp(-q*T[NT-1])*BS.NormCDF(g);

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
            double dPdxdy   = K*Math.Exp(-rf*T[NT-1])*PHIfxy   - Math.Exp(-q*T[NT-1])*S*(PHIgy + PHIgxy);
            double dPdx2dy  = K*Math.Exp(-rf*T[NT-1])*PHIfx2y  - Math.Exp(-q*T[NT-1])*S*(PHIgy + 2.0*PHIgxy + PHIgx2y);
            double dPdy2    = K*Math.Exp(-rf*T[NT-1])*PHIfy2   - Math.Exp(-q*T[NT-1])*S*PHIgy2;
            double dPdx2dy2 = K*Math.Exp(-rf*T[NT-1])*PHIfx2y2 - Math.Exp(-q*T[NT-1])*S*(PHIgy + 2.0*PHIgxy + PHIgx2y + PHIgy2 + 2.0*PHIgxy2 + PHIgx2y2);

            // Benhamou, Gobet, Miri expansion
            double Put = BSPut + A1T*dPdxdy + A2T*dPdx2dy + B0T*dPdy2 + B2T*dPdx2dy2;

            // Return the put or the call by put-call parity
            if(PutCall == "P")
                return Put;
            else
                return Put - K*Math.Exp(-rf*T[NT-1]) + S*Math.Exp(-q*T[NT-1]);
        }
    }
}


