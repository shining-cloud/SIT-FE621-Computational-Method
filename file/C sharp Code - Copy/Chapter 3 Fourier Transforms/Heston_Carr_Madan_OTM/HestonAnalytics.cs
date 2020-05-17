using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Carr_Madan_OTM
{
    class HestonPrice
    {
        // Heston Integrand
        public double HestonProb(double phi,double kappa,double theta,double lambda,double rho,double sigma,double T,
                          double K,double S,double r,double v0,int Pnum,int Trap)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            double x = Math.Log(S);
            double a = kappa * theta;
            Complex b,u,d,g,c,D,G,C,f,integrand = new Complex();

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
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2.0) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                D = (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*T))/(1.0-c*Complex.Exp(-d*T)));
                G = (1.0 - c*Complex.Exp(-d*T))/(1.0-c);
                C = r*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi - d)*T - 2.0*Complex.Log(G));
            }
            else
            {
                // Original Heston formulation.
                G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
                C = r*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi + d)*T - 2.0*Complex.Log(G));
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
                          double K,double S,double r,double v0,int Trap)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            double x = Math.Log(S);
            double a = kappa * theta;
            Complex b,u,d,g,c,D,G,C,f = new Complex();

            u = -0.5;
            b = kappa + lambda;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2.0) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                D = (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*T))/(1.0-c*Complex.Exp(-d*T)));
                G = (1.0 - c*Complex.Exp(-d*T))/(1.0-c);
                C = r*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi - d)*T - 2.0*Complex.Log(G));
            }
            else
            {
                // Original Heston formulation.
                G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
                C = r*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi + d)*T - 2.0*Complex.Log(G));
                D = (b - rho*sigma*i*phi + d)/sigma/sigma*((1.0-Complex.Exp(d*T))/(1.0-g*Complex.Exp(d*T)));
            }

            // The characteristic function.
            return Complex.Exp(C + D*v0 + i*phi*x);
        }

        // Heston Price by Gauss-Laguerre Integration
        public double HestonPriceGaussLaguerre(string Integrand,string PutCall,double S,double K,double r,double T,
                                               double kappa,double theta,double sigma,double v0,double lambda,double rho,
                                               double[] x,double[] w,int trap,double alpha)
        {
            double pi = Math.PI;
            if(Integrand == "Heston")
            {
                double[] int1 = new Double[32];
                double[] int2 = new Double[32];

                // Numerical integration
                for(int k=0;k<=31;k++)
                {
                    int1[k] = w[k] * HestonProb(x[k],kappa,theta,lambda,rho,sigma,T,K,S,r,v0,1,trap);
                    int2[k] = w[k] * HestonProb(x[k],kappa,theta,lambda,rho,sigma,T,K,S,r,v0,2,trap);
                }

                // Define P1 and P2
                double P1 = 0.5 + 1.0/pi*int1.Sum();
                double P2 = 0.5 + 1.0/pi*int2.Sum();

                // The call price
                double HestonC = S*P1 - K*Math.Exp(-r*T)*P2;

                // The put price by put-call parity
                double HestonP = HestonC - S + K*Math.Exp(-r*T);

                // Output the option price
                if(PutCall == "C")
                    return HestonC;
                else
                    return HestonP;
            }
            else if(Integrand == "CarrMadan")
            {
                double[] int1 = new Double[32];
                for(int k=0;k<=31;k++)
                    int1[k] = w[k] * CarrMadanIntegrandOTM(x[k],kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap);

                return int1.Sum() / pi;
            }
            else if(Integrand == "CarrMadanDamped")
            {
                double[] int1 = new Double[32];
                for(int k=0;k<=31;k++)
                    int1[k] = w[k] * CarrMadanDampedIntegrandOTM(x[k],kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,alpha);

                return 1.0/Math.Sinh(alpha*Math.Log(K))*int1.Sum() / pi;
            }
            else return 0.0;
        }

        // Returns the undamped Carr-Madan integrand for OTM options
        public double CarrMadanIntegrandOTM(double u,double kappa,double theta,double lambda,double rho,double sigma,
                                         double T,double K,double S,double r,double v0,int trap)
        {
            Complex i   = new Complex(0.0,1.0);
            Complex phi = HestonCF(u-i,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap);
            Complex integrand = Complex.Exp(-i*u*Math.Log(K)) * Math.Exp(-r*T)
                              * (Complex.Pow(S,(i*u+1.0))/(1.0+i*u) 
                              - Math.Exp(r*T)*Complex.Pow(S,(i*u+1.0))/(i*u) - phi/(u*u-i*u));
            return integrand.Real;
        }
        // Returns the damped Carr-Madan integrand for OTM options
        public double CarrMadanDampedIntegrandOTM(double v,double kappa,double theta,double lambda,double rho,double sigma,
                                                  double T,double K,double S,double r,double v0,int Trap,double alpha)
        {
            // Calculate z(v-ia)
            Complex i = new Complex(0.0,1.0);
            Complex u = v - i*alpha;
            Complex phi = HestonCF(u-i,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,Trap);
            Complex z1  = Math.Exp(-r*T)*(Complex.Pow(S,(i*u+1.0))/(1.0+i*u) - Math.Exp(r*T)*Complex.Pow(S,(i*u+1.0))/(i*u) - phi/(u*u-i*u));

            // Calculate z(v+ia)
            u  =  v + i*alpha;
            phi = HestonCF(u-i,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,Trap);
            Complex z2  = Math.Exp(-r*T)*(Complex.Pow(S,(i*u+1.0))/(1.0+i*u) - Math.Exp(r*T)*Complex.Pow(S,(i*u+1.0))/(i*u) - phi/(u*u-i*u));

            // Calculate the Fourier transform of y 
            Complex y = Complex.Exp(-i*u*Math.Log(K)) * (z1 - z2)/2.0;

            // Return the real part only
            return y.Real;
        }
    }
}

