using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Heston
{
    class MSExpansionBS
    {
        public double MSPutBS(double y,MSsetBS mssettings)
        {
            // The settings and parameters
            double theta = mssettings.theta;
            double K     = mssettings.K;
            double sigma = mssettings.sigma;
            double r     = mssettings.r;
            double q     = mssettings.q;
            double T     = mssettings.T;

            // The drift
            double mu = r-q;

            // The "C" coefficients evaluated at theta = y
            
            BisectionPricing BP = new BisectionPricing();
            double cdf = BP.NormCDF(y);
            double pdf = BP.NormPDF(y);
            double C1 = sigma*y*K/(y*cdf+pdf);
            double C2 = -1.0/2.0*(C1*cdf*Math.Pow(sigma,2.0)-2.0*C1*cdf*mu+Math.Pow(sigma,3.0)*Math.Pow(y,2.0)*K)/sigma/(cdf*Math.Pow(y,2.0)+cdf+y*pdf);
            double C3 =  1.0/24.0*(-24.0*y*cdf*Math.Pow(sigma,3.0)*C2+48.0*y*cdf*sigma*C2*mu+24.0*y*cdf*Math.Pow(sigma,2.0)*r*C1-24.0*pdf*C2*Math.Pow(sigma,3.0)+48.0*pdf*C2*mu*sigma+24.0*pdf*r*C1*Math.Pow(sigma,2.0)-3.0*pdf*C1*Math.Pow(sigma,4.0)+12.0*pdf*C1*Math.Pow(sigma,2.0)*mu-12.0*pdf*C1*Math.Pow(mu,2.0)+4.0*Math.Pow(sigma,5.0)*Math.Pow(y,3.0)*K)/Math.Pow(sigma,2.0)/(cdf*Math.Pow(y,3.0)+3.0*y*cdf+pdf*Math.Pow(y,2.0)+2.0*pdf);
            double C4 = -1.0/48.0*(-48.0*cdf*Math.Pow(sigma,3.0)*Math.Pow(y,2.0)*r*C2+72.0*cdf*Math.Pow(sigma,4.0)*Math.Pow(y,2.0)*C3-144.0*cdf*Math.Pow(sigma,2.0)*Math.Pow(y,2.0)*C3*mu-48.0*cdf*Math.Pow(sigma,3.0)*r*C2+72.0*cdf*Math.Pow(sigma,4.0)*C3-144.0*cdf*Math.Pow(sigma,2.0)*C3*mu+12.0*cdf*Math.Pow(sigma,5.0)*C2-48.0*cdf*Math.Pow(sigma,3.0)*C2*mu-24.0*cdf*Math.Pow(sigma,4.0)*r*C1+48.0*cdf*sigma*C2*Math.Pow(mu,2.0)+48.0*cdf*Math.Pow(sigma,2.0)*mu*r*C1-48.0*y*pdf*r*C2*Math.Pow(sigma,3.0)+72.0*y*pdf*C3*Math.Pow(sigma,4.0)-144.0*y*pdf*C3*mu*Math.Pow(sigma,2.0)-y*pdf*Math.Pow(sigma,6.0)*C1+6.0*y*pdf*Math.Pow(sigma,4.0)*C1*mu-12.0*y*pdf*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,2.0)+8.0*y*pdf*C1*Math.Pow(mu,3.0)+2.0*Math.Pow(sigma,7.0)*Math.Pow(y,4.0)*K)/Math.Pow(sigma,3.0)/(cdf*Math.Pow(y,4.0)+6.0*cdf*Math.Pow(y,2.0)+3.0*cdf+pdf*Math.Pow(y,3.0)+5.0*y*pdf);
            double C5 =  1.0/1920.0*(16.0*Math.Pow(sigma,9.0)*Math.Pow(y,5.0)*K-80.0*pdf*Math.Pow(sigma,7.0)*C2+5.0*pdf*Math.Pow(sigma,8.0)*C1+80.0*pdf*C1*Math.Pow(mu,4.0)-640.0*pdf*Math.Pow(r,2.0)*C3*Math.Pow(sigma,4.0)*C1+320.0*pdf*r*C3*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,2.0)+640.0*pdf*Math.Pow(mu,2.0)*r*C1*Math.Pow(sigma,2.0)+640.0*pdf*r*C3*Math.Pow(sigma,5.0)*C2+640.0*pdf*sigma*C2*Math.Pow(mu,3.0)+480.0*pdf*Math.Pow(sigma,5.0)*C2*mu-960.0*pdf*Math.Pow(sigma,3.0)*C2*Math.Pow(mu,2.0)+160.0*pdf*Math.Pow(sigma,6.0)*r*C1+120.0*pdf*Math.Pow(sigma,4.0)*C1*Math.Pow(mu,2.0)-40.0*pdf*Math.Pow(sigma,6.0)*C1*mu-160.0*pdf*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,3.0)+3840.0*pdf*r*Math.Pow(sigma,4.0)*C3+1280.0*pdf*r*Math.Pow(sigma,5.0)*C2-320.0*pdf*Math.Pow(r,2.0)*Math.Pow(sigma,4.0)*C1+15360.0*pdf*Math.Pow(sigma,3.0)*C4*mu+5760.0*pdf*Math.Pow(sigma,4.0)*C3*mu-5760.0*pdf*Math.Pow(sigma,2.0)*C3*Math.Pow(mu,2.0)+23040.0*y*cdf*Math.Pow(sigma,3.0)*C4*mu+1920.0*y*cdf*Math.Pow(sigma,5.0)*r*C2-3840.0*y*cdf*Math.Pow(sigma,3.0)*r*C2*mu-3840.0*pdf*Math.Pow(y,2.0)*Math.Pow(sigma,5.0)*C4-5.0*pdf*Math.Pow(y,2.0)*Math.Pow(sigma,8.0)*C1-80.0*pdf*Math.Pow(y,2.0)*C1*Math.Pow(mu,4.0)-960.0*y*cdf*Math.Pow(sigma,4.0)*Math.Pow(r,2.0)*C1+5760.0*y*cdf*Math.Pow(sigma,4.0)*C3*mu-5760.0*y*cdf*Math.Pow(sigma,2.0)*C3*Math.Pow(mu,2.0)+40.0*pdf*Math.Pow(y,2.0)*Math.Pow(sigma,6.0)*C1*mu-120.0*pdf*Math.Pow(y,2.0)*Math.Pow(sigma,4.0)*C1*Math.Pow(mu,2.0)+160.0*pdf*Math.Pow(y,2.0)*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,3.0)+7680.0*pdf*Math.Pow(y,2.0)*Math.Pow(sigma,3.0)*C4*mu+1920.0*pdf*Math.Pow(y,2.0)*r*Math.Pow(sigma,4.0)*C3-3840.0*Math.Pow(y,3.0)*cdf*Math.Pow(sigma,5.0)*C4-11520.0*y*cdf*Math.Pow(sigma,5.0)*C4-1440.0*y*cdf*Math.Pow(sigma,6.0)*C3-320.0*pdf*r*C3*Math.Pow(sigma,4.0)*C1*mu-1280.0*pdf*r*C3*Math.Pow(sigma,3.0)*C2*mu-640.0*pdf*Math.Pow(sigma,4.0)*mu*r*C1+80.0*pdf*r*C3*Math.Pow(sigma,6.0)*C1-7680.0*pdf*Math.Pow(sigma,5.0)*C4-1440.0*pdf*Math.Pow(sigma,6.0)*C3+1920.0*Math.Pow(y,3.0)*cdf*Math.Pow(sigma,4.0)*r*C3+7680.0*Math.Pow(y,3.0)*cdf*Math.Pow(sigma,3.0)*C4*mu+5760.0*y*cdf*Math.Pow(sigma,4.0)*r*C3-2560.0*pdf*r*Math.Pow(sigma,3.0)*C2*mu)/Math.Pow(sigma,4.0)/(cdf*Math.Pow(y,5.0)+10.0*cdf*Math.Pow(y,3.0)+15.0*y*cdf+pdf*Math.Pow(y,4.0)+9.0*pdf*Math.Pow(y,2.0)+8.0*pdf);

            //Set 1 polynomials
            double P01 = theta;
            double P11 = 0.0;
            double Q01 = 1.0;
            double Q11 = 0.0;

            //Set 2 polynomials
            double P02 = Math.Pow(theta,2.0)+1.0;
            double P12 = -1.0/2.0*C1*(-Math.Pow(sigma,2.0)+2.0*mu)/sigma;
            double Q02 = theta;
            double Q12 = 0.0;

            //Set 3 polynomials
            double P03 = Math.Pow(theta,3.0)+3.0*theta;
            double P13 = -theta*(-C2*Math.Pow(sigma,2.0)+2.0*C2*mu+r*C1*sigma)/sigma;
            double Q03 = Math.Pow(theta,2.0)+2.0;
            double Q13 = -1.0/8.0*(-8.0*C2*Math.Pow(sigma,3.0)+16.0*C2*mu*sigma+8.0*r*C1*Math.Pow(sigma,2.0)-C1*Math.Pow(sigma,4.0)+4.0*C1*Math.Pow(sigma,2.0)*mu-4.0*C1*Math.Pow(mu,2.0))/Math.Pow(sigma,2.0);

            //Set 4 polynomials
            double P04 = Math.Pow(theta,4.0)+6.0*Math.Pow(theta,2.0)+3.0;
            double P14 = -1.0/2.0*(2.0*r*C2*sigma-3.0*C3*Math.Pow(sigma,2.0)+6.0*C3*mu)/sigma*Math.Pow(theta,2.0)+1.0/4.0*(-4.0*r*C2*Math.Pow(sigma,2.0)+6.0*C3*Math.Pow(sigma,3.0)-12.0*C3*sigma*mu+Math.Pow(sigma,4.0)*C2-4.0*Math.Pow(sigma,2.0)*C2*mu-2.0*Math.Pow(sigma,3.0)*r*C1+4.0*C2*Math.Pow(mu,2.0)+4.0*mu*r*C1*sigma)/Math.Pow(sigma,2.0);
            double Q04 = Math.Pow(theta,3.0)+5.0*theta;
            double Q14 = -1.0/48.0*(-72.0*C3*Math.Pow(sigma,4.0)+144.0*C3*mu*Math.Pow(sigma,2.0)+48.0*r*C2*Math.Pow(sigma,3.0)+Math.Pow(sigma,6.0)*C1-6.0*Math.Pow(sigma,4.0)*C1*mu+12.0*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,2.0)-8.0*C1*Math.Pow(mu,3.0))/Math.Pow(sigma,3.0)*theta;

            //Set 5 polynomials
            double P05 = Math.Pow(theta,5.0)+10.0*Math.Pow(theta,3.0)+15.0*theta;
            double P15 = -(r*C3*sigma-2.0*C4*Math.Pow(sigma,2.0)+4.0*C4*mu)/sigma*Math.Pow(theta,3.0)+1.0/4.0*(-12.0*r*C3*Math.Pow(sigma,2.0)+24.0*C4*Math.Pow(sigma,3.0)-48.0*C4*sigma*mu-4.0*r*C2*Math.Pow(sigma,3.0)+8.0*r*sigma*C2*mu+2.0*Math.Pow(r,2.0)*Math.Pow(sigma,2.0)*C1+3.0*C3*Math.Pow(sigma,4.0)-12.0*C3*mu*Math.Pow(sigma,2.0)+12.0*C3*Math.Pow(mu,2.0))/Math.Pow(sigma,2.0)*theta;
            double Q05 = Math.Pow(theta,4.0)+9.0*Math.Pow(theta,2.0)+8;
            double Q15 = -1.0/384.0*(-Math.Pow(sigma,8.0)*C1+8.0*Math.Pow(sigma,6.0)*C1*mu-24.0*Math.Pow(sigma,4.0)*C1*Math.Pow(mu,2.0)+32.0*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,3.0)-16.0*C1*Math.Pow(mu,4.0)-768.0*C4*Math.Pow(sigma,5.0)+1536.0*C4*Math.Pow(sigma,3.0)*mu+384.0*r*C3*Math.Pow(sigma,4.0))/Math.Pow(sigma,4.0)*Math.Pow(theta,2.0)+1.0/384.0*(-128.0*r*C3*Math.Pow(sigma,5.0)*C2-128.0*Math.Pow(mu,2.0)*r*C1*Math.Pow(sigma,2.0)+1152.0*C3*Math.Pow(mu,2.0)*Math.Pow(sigma,2.0)+32.0*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,3.0)-24.0*Math.Pow(sigma,4.0)*C1*Math.Pow(mu,2.0)+8.0*Math.Pow(sigma,6.0)*C1*mu-256.0*Math.Pow(sigma,5.0)*r*C2-1152.0*Math.Pow(sigma,4.0)*C3*mu-768.0*r*C3*Math.Pow(sigma,4.0)+128.0*Math.Pow(r,2.0)*C3*Math.Pow(sigma,4.0)*C1-3072.0*C4*Math.Pow(sigma,3.0)*mu+256.0*r*C3*Math.Pow(sigma,3.0)*C2*mu-16.0*r*C3*Math.Pow(sigma,6.0)*C1+64.0*r*C3*Math.Pow(sigma,4.0)*C1*mu-64.0*r*C3*Math.Pow(sigma,2.0)*C1*Math.Pow(mu,2.0)+64.0*Math.Pow(sigma,4.0)*Math.Pow(r,2.0)*C1-128.0*sigma*C2*Math.Pow(mu,3.0)+192.0*Math.Pow(sigma,3.0)*C2*Math.Pow(mu,2.0)-32.0*Math.Pow(sigma,6.0)*r*C1-96.0*Math.Pow(sigma,5.0)*C2*mu+1536.0*C4*Math.Pow(sigma,5.0)+512.0*mu*r*C2*Math.Pow(sigma,3.0)-16.0*C1*Math.Pow(mu,4.0)+288.0*Math.Pow(sigma,6.0)*C3-Math.Pow(sigma,8.0)*C1+16.0*Math.Pow(sigma,7.0)*C2+128.0*Math.Pow(sigma,4.0)*mu*r*C1)/Math.Pow(sigma,4.0);

            // The Black-Scholes American put approximation
            cdf = BP.NormCDF(theta);
            pdf = BP.NormPDF(theta);

            double Price = (C1*(P01*cdf + Q01*pdf) + P11*cdf + Q11*pdf)*Math.Pow(T,0.5)
                  + (C2*(P02*cdf + Q02*pdf) + P12*cdf + Q12*pdf)*T
                  + (C3*(P03*cdf + Q03*pdf) + P13*cdf + Q13*pdf)*Math.Pow(T,1.5)
                  + (C4*(P04*cdf + Q04*pdf) + P14*cdf + Q14*pdf)*Math.Pow(T,2.0)
                  + (C5*(P05*cdf + Q05*pdf) + P15*cdf + Q15*pdf)*Math.Pow(T,2.5);

            return Price;
        }
    }
}
