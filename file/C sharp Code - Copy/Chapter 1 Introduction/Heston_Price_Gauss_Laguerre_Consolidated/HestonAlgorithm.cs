using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Price_Gauss_Laguerre_Consolidated
{
    class HestonPriceConsolidated
    {
        // Heston Integrand
        public double HestonProbConsol(double phi,HParam param,OpSet settings)
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
            Complex b1,u1,d1,g1,c1,D1,G1,C1,f1,
                    b2,u2,d2,g2,c2,D2,G2,C2,f2,integrand = new Complex();

            // The first characteristic function
            u1 = 0.5;
            b1 = kappa + lambda - rho*sigma;
            d1 = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b1,2) - sigma*sigma*(2.0*u1*i*phi - phi*phi));
            g1 = (b1 - rho*sigma*i*phi + d1) / (b1 - rho*sigma*i*phi - d1);
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c1 = 1.0/g1;
                D1 = (b1 - rho*sigma*i*phi - d1)/sigma/sigma*((1.0-Complex.Exp(-d1*T))/(1.0-c1*Complex.Exp(-d1*T)));
                G1 = (1.0 - c1*Complex.Exp(-d1*T))/(1.0-c1);
                C1 = (r-q)*i*phi*T + a/sigma/sigma*((b1 - rho*sigma*i*phi - d1)*T - 2.0*Complex.Log(G1));
            }
            else
            {
                // Original Heston formulation.
                G1 = (1.0 - g1*Complex.Exp(d1*T))/(1.0-g1);
                C1 = (r-q)*i*phi*T + a/sigma/sigma*((b1 - rho*sigma*i*phi + d1)*T - 2.0*Complex.Log(G1));
                D1 = (b1 - rho*sigma*i*phi + d1)/sigma/sigma*((1.0-Complex.Exp(d1*T))/(1.0-g1*Complex.Exp(d1*T)));
            }
            f1 = Complex.Exp(C1 + D1*v0 + i*phi*x);

            // The second characteristic function
            u2 = -0.5;
            b2 = kappa + lambda;
            d2 = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b2,2) - sigma*sigma*(2.0*u2*i*phi - phi*phi));
            g2 = (b2 - rho*sigma*i*phi + d2) / (b2 - rho*sigma*i*phi - d2);
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c2 = 1.0/g2;
                D2 = (b2 - rho*sigma*i*phi - d2)/sigma/sigma*((1.0-Complex.Exp(-d2*T))/(1.0-c2*Complex.Exp(-d2*T)));
                G2 = (1.0 - c2*Complex.Exp(-d2*T))/(1.0-c2);
                C2 = (r-q)*i*phi*T + a/sigma/sigma*((b2 - rho*sigma*i*phi - d2)*T - 2.0*Complex.Log(G2));
            }
            else
            {
                // Original Heston formulation.
                G2 = (1.0 - g2*Complex.Exp(d2*T))/(1.0-g2);
                C2 = (r-q)*i*phi*T + a/sigma/sigma*((b2 - rho*sigma*i*phi + d2)*T - 2.0*Complex.Log(G2));
                D2 = (b2 - rho*sigma*i*phi + d2)/sigma/sigma*((1.0-Complex.Exp(d2*T))/(1.0-g2*Complex.Exp(d2*T)));
            }
            f2 = Complex.Exp(C2 + D2*v0 + i*phi*x);

            // The integrand.
            integrand = Complex.Exp(-i*phi*Complex.Log(K))/i/phi*(S*Complex.Exp(-q*T)*f1 - K*Complex.Exp(-r*T)*f2);

            // Return the real part of the integrand.
            return integrand.Real;
        }

        // Heston Price by Gauss-Laguerre Integration
        public double HestonPriceConsol(HParam param,OpSet settings,double[] x,double[] w)
        {
            double[] int1 = new Double[32];
            // Numerical integration
            for(int j=0;j<=31;j++)
            {
                int1[j] = w[j] * HestonProbConsol(x[j],param,settings);
            }

            // Define P1 and P2
            double pi = Math.PI;
            double I = int1.Sum();

            // The call price
            double S = settings.S;
            double K = settings.K;
            double r = settings.r;
            double q = settings.q;
            double T = settings.T;
            string PutCall = settings.PutCall;
            double HestonC = 0.5*S*Math.Exp(-q*T) - 0.5*K*Math.Exp(-r*T) + I/pi;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);

            // Output the option price
            if(PutCall == "C")
                return HestonC;
            else
                return HestonP;

        }
    }
}

