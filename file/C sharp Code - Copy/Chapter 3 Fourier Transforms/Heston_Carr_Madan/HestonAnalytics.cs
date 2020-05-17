using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Carr_Madan
{
    class HestonPrice
    {
        // Heston Integrand
        public double HestonProb(double phi,double kappa,double theta,double lambda,double rho,double sigma,double T,
                          double K,double S,double r,double q,double v0,int Pnum,int Trap)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            Complex b,u,d,g,c,D,G,C,f,integrand = new Complex();
            double x = Math.Log(S);
            double a = kappa * theta;

            // Parameters "u" and "b" are different for P1 and P2
            if(Pnum==1)
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
            if(Trap==1)
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

            // The real part of the integrand
            return integrand.Real;
        }

        // Heston Characteristic Function (f2)
        public Complex HestonCF(Complex phi,double kappa,double theta,double lambda,double rho,double sigma,double T,
                          double K,double S,double r,double q,double v0,int Trap)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            Complex b,u,d,g,c,D,G,C = new Complex();
            double x = Math.Log(S);
            double a = kappa * theta;

            u = -0.5;
            b = kappa + lambda;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(Trap==1)
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
            return Complex.Exp(C + D*v0 + i*phi*x);
        }

        // Heston Price by Gauss-Laguerre Integration
        public double HestonPriceGaussLaguerre(string Integrand,string PutCall,double alpha,double S,double K,double r,double q,double T,
                                               double kappa,double theta,double sigma,double v0,double lambda,double rho,double[] x,double[] w,int trap)
        {
            if(Integrand == "Heston")
            {
                double[] int1 = new Double[32];
                double[] int2 = new Double[32];

                // Numerical integration
                for(int k=0;k<=31;k++)
                {
                    int1[k] = w[k] * HestonProb(x[k],kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
                    int2[k] = w[k] * HestonProb(x[k],kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
                }

                // Define P1 and P2
                double pi = Math.PI;
                double P1 = 0.5 + 1.0/pi*int1.Sum();
                double P2 = 0.5 + 1.0/pi*int2.Sum();

                // The call price
                double HestonC = S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2;

                // The put price by put-call parity
                double HestonP = HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);

                // Output the option price
                if(PutCall == "C")
                    return HestonC;
                else
                    return HestonP;
            }
            else
            {
                double[] int1 = new Double[32];
                for(int k=0;k<=31;k++)
                    int1[k] = w[k] * CarrMadanIntegrand(x[k],alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,PutCall);

                // The Option Price
                double pi = Math.PI;
                if(PutCall == "C")
                    return Math.Exp(-alpha*Math.Log(K)) * int1.Sum() / pi;
                else
                    return Math.Exp(alpha*Math.Log(K)) * int1.Sum() / pi;
                return int1.Sum();
            }
        }

        // Returns the Carr-Madan integrand,
        // based on the Heston characteristic function, f2 (the second CF)
        public double CarrMadanIntegrand(double u,double alpha,double kappa,double theta,double lambda,double rho,double sigma,
                                         double T,double K,double S,double r,double q,double v0,int trap,string PutCall)
        {
            Complex i   = new Complex(0.0,1.0);       // Imaginary unit
            Complex one = new Complex(1.0,0.0);
            Complex two = new Complex(2.0,0.0);
            Complex integrand;
            if(PutCall == "C")
                integrand = Complex.Exp(-i*u*Complex.Log(K)) * Complex.Exp(-r*T)
		                  * HestonCF(u-(alpha+one)*i,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap)
		                  / (alpha*alpha + alpha - u*u + i*(two*alpha+one)*u);
            else
                integrand = Complex.Exp(-i*u*Complex.Log(K)) * Complex.Exp(-r*T)
		                  * HestonCF(u-(-alpha+one)*i,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap)
                 		  / (alpha*alpha - alpha - u*u + i*(-two*alpha+one)*u);
            return integrand.Real;
        }
    }
}

