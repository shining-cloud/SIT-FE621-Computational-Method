using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_2001_Greeks
{
    class LewisGreeks
    {
        // Integrand for Greeks from Equation 3.11 of Lewis (2001)
        public double LewisIntegrand311Greeks(Complex u,double kappa,double theta,double lambda,double rho,double sigma,double tau,double S,double K,double r,double q,double v0,int Trap,string Greek)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            // Change u to u - i/2
            u -= i/2.0;

            // Log of the stock price.
            double x0 = Math.Log(S);

            // Heston a and b parameters
            double a = kappa*theta;
            double b = kappa + lambda;
            Complex d = Complex.Sqrt(Complex.Pow(rho*sigma*i*u - b,2.0) + sigma*sigma*(i*u + u*u));
            Complex g = (b - rho*sigma*i*u + d) / (b - rho*sigma*i*u - d);
            Complex c=0.0,D=0.0,G=0.0,C=0.0;
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                D = (b - rho*sigma*i*u - d)/sigma/sigma*((1.0-Complex.Exp(-d*tau))/(1.0-c*Complex.Exp(-d*tau)));
                G = (1.0 - c*Complex.Exp(-d*tau))/(1.0-c);
                C = (r-q)*i*u*tau + a/sigma/sigma*((b - rho*sigma*i*u - d)*tau - 2.0*Complex.Log(G));
            }
            else if(Trap==0)
            {
                // Original Heston formulation.
                G = (1.0 - g*Complex.Exp(d*tau))/(1.0-g);
                C = (r-q)*i*u*tau + a/sigma/sigma*((b - rho*sigma*i*u + d)*tau - 2.0*Complex.Log(G));
                D = (b - rho*sigma*i*u + d)/sigma/sigma*((1.0-Complex.Exp(d*tau))/(1.0-g*Complex.Exp(d*tau)));
            }

            // The Heston characteristic function f2 for ln S(T)
            Complex f2 = Complex.Exp(C + D*v0 + i*u*x0);

            // Derivative of f2 w.r.t. tau
            Complex dD = d*Complex.Exp(d*tau)*(b-rho*sigma*u*i+d)*(g-1.0)/sigma/sigma/Complex.Pow(1-g*Complex.Exp(d*tau),2.0);
            Complex dC = (r-q)*u*i + kappa*theta/sigma/sigma * ((b-rho*sigma*u*i+d) + 2.0*g*d*Complex.Exp(d*tau)/(1.0-g*Complex.Exp(d*tau)));
            Complex dt = (dC + dD*v0);

            // Derivative of f2 w.r.t. r
            Complex dr = i*u*tau;

            // Change u - i/2 back to u
            u += i/2;

            // The integrand
            Complex y = 0.0;
            if(Greek == "Price")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2;
            else if(Greek == "Delta")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2*(i*u+0.5)/S;
            else if(Greek == "Gamma")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2*(i*u+0.5)*(i*u-0.5)/S/S;
            else if(Greek == "Rho")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2*(-tau+dr);
            else if(Greek == "Theta")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2*(-r+dt);
            else if(Greek == "Vega1")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2*D;
            else if(Greek == "Volga")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2*(4.0*D*D*v0 + 2.0*D);
            else if(Greek == "Vanna")
                y = Complex.Pow(K,-i*u)/(u*u+0.25)*Math.Exp(-r*tau)*f2*D*(i*u+0.5)/S;
            return y.Real;
        }
        // Return the Greeks
        public double LewisGreeks311(double S,double K,double r,double q,double T,double theta,double kappa,double sigma,double rho,double v0,int trap,double[] x,double[] w,string Greek)
        {
            double lambda = 0.0;

            // Compute the integral.
            // The integration variable is complex k = kr + (1/2)*i
            int Nx = x.Length;
            Complex u;
            double[] Int = new double[Nx];
            for(int j=0;j<=Nx-1;j++)
            {
                u = x[j];
                Int[j] = w[j] * LewisIntegrand311Greeks(u,kappa,theta,lambda,rho,sigma,T,S,K,r,q,v0,trap,Greek);
            }
            double Integral = Int.Sum();
            double pi = Math.PI;

            // Equation (3.11) in Lewis (2011)
            if(Greek == "Price")
                return S*Math.Exp(-q*T) - Math.Sqrt(K)/pi*Integral;
            else if(Greek == "Delta")
                return Math.Exp(-q*T) - Math.Sqrt(K)/pi*Integral;
            else if((Greek == "Gamma") || (Greek == "Rho"))
                return -Math.Sqrt(K)/pi*Integral;
            else if(Greek == "Theta")
                return q*S*Math.Exp(-q*T) + Math.Sqrt(K)/pi*Integral;
            else if((Greek == "Vega1") || (Greek == "Vanna"))
                return -Math.Sqrt(K)/pi*Integral*2*Math.Sqrt(v0);
            else if(Greek == "Volga")
                return -Math.Sqrt(K)/pi*Integral;
            else
                return 0.0;
        }
    }
}
