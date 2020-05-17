using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Euler_and_Milstein
{
    class HestonPrice
    {
        // Heston Integrand
        public double HestonProb(double phi,HParam param, OpSet settings, int Pnum)
        {
            Complex i  = new Complex(0.0,1.0);                   // Imaginary unit
            double S = settings.S;
            double K = settings.K;
            double T = settings.T;
            double r = settings.r;
            double q = settings.q;
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0 = param.v0;
            double rho = param.rho;
            double lambda = param.lambda;
            double x = Math.Log(S);
            double a = kappa*theta;
            int Trap = settings.trap;
            Complex b,u,d,g,c,D,G,C,f,integrand = new Complex();

            // Parameters "u" and "b" are different for P1 and P2
            if (Pnum==1)
            {
            	u = 0.5;
	            b = kappa + lambda - rho*sigma;
            }
            else
            {
	            u = -0.5;
	            b = kappa + lambda;
            }
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if (Trap==1)
            {
            	// "Little Heston Trap" formulation
    	        c = 1.0/g;
        	    D = (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*T))/(1.0-c*Complex.Exp(-d*T)));
    	        G = (1.0 - c*Complex.Exp(-d*T))/(1-c);
        	    C = (r-q)*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi - d)*T - 2.0*Complex.Log(G));
            }
            else
            {
                // Original Heston formulation.
        	    G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
            	C = (r-q)*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi + d)*T - 2.0*Complex.Log(G));
	            D = (b - rho*sigma*i*phi + d)/sigma/sigma*((1.0-Complex.Exp(d*T))/(1.0-g*Complex.Exp(d*T)));
            }

            // The characteristic function.
            f = Complex.Exp(C + D*v0 + i*phi*x);

            // The integrand.
            integrand = Complex.Exp(-i*phi*Complex.Log(K))*f/i/phi;

            // Return the real part of the integrand.
            return integrand.Real;
        }

        // Heston Price by Gauss-Laguerre Integration
        public double HestonPriceGaussLaguerre(HParam param, OpSet settings,double[] x, double[] w)
        {
            double[] int1 = new Double[32];
            double[] int2 = new Double[32];
            // Numerical integration
            for (int j=0; j<=31; j++)
            {
                int1[j] = w[j] * HestonProb(x[j],param,settings,1);
	            int2[j] = w[j] * HestonProb(x[j],param,settings,2);
            }

            // Define P1 and P2
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*int1.Sum();
            double P2 = 0.5 + 1.0/pi*int2.Sum();

            // The call price
            double S = settings.S;
            double K = settings.K;
            double T = settings.T;
            double r = settings.r;
            double q = settings.q;
            string PutCall = settings.PutCall;
            double HestonC = S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);

            // Output the option price
            if (PutCall == "C")
	            return HestonC;
            else
	            return HestonP;
        }
    }
}

