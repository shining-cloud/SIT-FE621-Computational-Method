using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Kahl_Jackel_Integration
{
    class KahlJackel
    {
        // Kahl and Jackel (2005) Heston price, using an integration range of [0,1].
        // From "Not-so-Complex Algorithms in the Heston model"
        public double HestonPriceKahlJackel(HParam param,OpSet settings,double[] X,double[] W)
        {
            int N = X.Length;
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            double S = settings.S;
            double K = settings.K;
            double T = settings.T;
            double r = settings.r;
            double q = settings.q;
            string PutCall = settings.PutCall;

            // Forward price
            double F = S*Math.Exp((r-q)*T);

            // Required quantities for integrand
            double Cinf = Math.Sqrt(1-rho*rho)/sigma*(v0 + kappa*theta*T);

            // Required quantities for f1
            double ImC1 = (Math.Exp((rho*sigma-kappa)*T)*theta*kappa + theta*kappa*((kappa-rho*sigma)*T-1.0))/2.0/Math.Pow(kappa-rho*sigma,2);
            double ImD1 = (1.0-Math.Exp((rho*sigma-kappa)*T))/2.0/(kappa-rho*sigma);

            // Required quantities for f2
            double ImC2 = -(Math.Exp(-kappa*T)*theta*kappa + theta*kappa*(kappa*T-1.0))/2.0/kappa/kappa;
            double ImD2 = -(1.0 - Math.Exp(-kappa*T/2.0))/2.0/kappa;

            double pi = Math.PI;
            double[] y = new double[N];
            double[] z = new double[N];
            HestonPrice HP = new HestonPrice();

            for(int u=0;u<=N-1;u++)
            {
                // Transformation of the abscissa from [-1,1] to [0,1]
                double x = 0.5*X[u] + 0.5;
                if(x == 0.0)
                    // Integrand at left abscissa 0
                    y[u] = 0.5*(F-K);
                else if(x == 1.0)
                {
                    // Integrand at right abscissa 1
                    double f1 = Math.Log(F/K) + ImC1 + ImD1*v0;
                    double f2 = Math.Log(F/K) + ImC2 + ImD2*v0;
                    y[u] = 0.5*(F-K) + (F*f1 - K*f1)/(pi*Cinf);
                }
                else
                {
                    // Integrand at remaining abscissas
                    double f1 = HP.HestonProb(-Math.Log(x)/Cinf,param,settings,1);
                    double f2 = HP.HestonProb(-Math.Log(x)/Cinf,param,settings,2);
                    y[u] = 0.5*(F-K) + (F*f1 - K*f2)/(x*pi*Cinf);
                }
                // Multiply by the weights
                z[u] = W[u] * y[u];
            }
            // The Call price
            double HCall = Math.Exp(-r*T)*0.5*z.Sum();

            // The put price by put call parity
            if(PutCall == "C")
                return HCall;
            else
                return HCall - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);
        }
    }
}
                /*
        // Kahl and Jackel (2005) Heston integrand
        // From "Not-so-Complex Algorithms in the Heston model"
        static double HestonIntegrandKahlJackel(double u,HParam param,OpSet settings)
        {
            double S  = settings.S;
            double K  = settings.K;
            double rf = settings.r;
            double q  = settings.q;
            double T  = settings.T;

            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double rho   = param.rho;
            double v0    = param.v0;

            // Forward price
            double F = S*Math.Exp((rf-q)*T);

            // Required quantities for integrand
            double Cinf = Math.Sqrt(1-rho*rho)/sigma*(v0 + kappa*theta*T);

            // Required quantities for f1
            double ImC1 = (Math.Exp((rho*sigma-kappa)*T)*theta*kappa + theta*kappa*((kappa-rho*sigma)*T-1.0))/2.0/Math.Pow(kappa-rho*sigma,2);
            double ImD1 = (1.0-Math.Exp((rho*sigma-kappa)*T))/2.0/(kappa-rho*sigma);

            // Required quantities for f2
            double ImC2 = -(Math.Exp(-kappa*T)*theta*kappa + theta*kappa*(kappa*T-1.0))/2.0/kappa/kappa;
            double ImD2 = - (1.0 - Math.Exp(-kappa*T/2.0))/2.0/kappa;

            double y;
            double pi = Math.PI;
            if (u == 0.0)
            	// Integrand at left abscissa 0
	            y = 0.5*(F-K);
            else if (u == 1.0)
            {
            	// Integrand at right abscissa 1
	            double f1 = Math.Log(F/K) + ImC1 + ImD1*v0;
	            double f2 = Math.Log(F/K) + ImC2 + ImD2*v0;
	            y = 0.5*(F-K) + (F*f1 - K*f1)/(pi*Cinf);
            }
            else
                {
            	// Integrand at remaining abscissas
	            double f1 = HestonProb(-Math.Log(u)/Cinf,param,settings,1);
            	double f2 = HestonProb(-Math.Log(u)/Cinf,param,settings,2);
            	y = 0.5*(F-K) + (F*f1 - K*f2)/(u*pi*Cinf);
            }
            return y;
        } */





