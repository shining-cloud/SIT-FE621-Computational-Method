using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Attari_Delta_and_Gamma
{
    class AttariPrice
    {
        // Attari Integrand
        public double AttariProb(double u,double kappa,double theta,double lambda,double rho,double sigma,double T,
                          double K,double S,double r,double q,double v0, int Trap)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            Complex b,d,g,c,D,G,C,f = new Complex();

            double a = kappa*theta;
            b = kappa + lambda;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*u - b,2.0) + sigma*sigma*(i*u + u*u));
            g = (b - rho*sigma*i*u + d) / (b - rho*sigma*i*u - d);
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                D = (b - rho*sigma*i*u - d)/sigma/sigma*((1.0-Complex.Exp(-d*T))/(1.0-c*Complex.Exp(-d*T)));
                G = (1.0 - c*Complex.Exp(-d*T))/(1-c);
                C = (r-q)*i*u*T + a/sigma/sigma*((b - rho*sigma*i*u - d)*T - 2.0*Complex.Log(G));
            }
            else
            {
                // Original Heston formulation.
                G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
                C = (r-q)*i*u*T + a/sigma/sigma*((b - rho*sigma*i*u + d)*T - 2.0*Complex.Log(G));
                D = (b - rho*sigma*i*u + d)/sigma/sigma*((1.0-Complex.Exp(d*T))/(1.0-g*Complex.Exp(d*T)));
            }

            // The Attari characteristic function.
            f = Complex.Exp(C + D*v0 - i*u*r*T);

            // The Attari (2004) integrand
            double L = Math.Log(Math.Exp(-r*T)*K/S);
            double integrand = ((f.Real + f.Imaginary/u)*Math.Cos(L*u) + (f.Imaginary - f.Real/u)*Math.Sin(L*u)) / (1 + u*u);
            return integrand;
        }

        // The Attari call price using Gauss-Laguerre 32 point integration+
        public double AttariPriceGaussLaguerre(string PutCall,double S,double K,double T,double r,double q,double kappa,double theta,double sigma,
                                               double lambda,double v0,double rho,int trap,double[] x,double[] w)
        {
            int N = x.Length;
            double[] int1 = new double[N];
            for(int k=0;k<=N-1;k++)
                int1[k] = w[k] * AttariProb(x[k],kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap);

            // The call price
            double pi = Math.PI;
            double HestonC = S*Math.Exp(-q*T) - K*Math.Exp(-r*T)*(0.5 + 1.0/pi*int1.Sum());

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


