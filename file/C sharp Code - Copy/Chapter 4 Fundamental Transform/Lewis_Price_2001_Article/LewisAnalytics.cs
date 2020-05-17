using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Price_2001_Article
{
    class LewisPrice
    {
        // Integrand in Equation 3.11 of Lewis (2001)
        public double LewisIntegrand311(Complex u,double kappa,double theta,double lambda,double rho,double sigma,double tau,double S,double K,double r,double q,double v0,int Trap)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
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
                D = (b - rho*sigma*i*u - d)/sigma/sigma*((1-Complex.Exp(-d*tau))/(1.0-c*Complex.Exp(-d*tau)));
                G = (1.0 - c*Complex.Exp(-d*tau))/(1.0-c);
                C = (r-q)*i*u*tau + a/sigma/sigma*((b - rho*sigma*i*u - d)*tau - 2.0*Complex.Log(G));
            }
            else if(Trap==0)
            {
                // Original Heston formulation.
                G = (1 - g*Complex.Exp(d*tau))/(1.0-g);
                C = (r-q)*i*u*tau + a/sigma/sigma*((b - rho*sigma*i*u + d)*tau - 2.0*Complex.Log(G));
                D = (b - rho*sigma*i*u + d)/sigma/sigma*((1.0-Complex.Exp(d*tau))/(1.0-g*Complex.Exp(d*tau)));
            }
            // The Heston characteristic function f2 for ln S(T)
            Complex CFlnST = Complex.Exp(C + D*v0 + i*u*x0);

            // The cf for the Levy process XT
            double Y = Math.Log(S) + (r-q)*tau;
            Complex CFXT = Complex.Exp(-i*u*Y)*CFlnST;

            // The integrand in Equation (3.11) of Lewis (2001)
            u += i/2;
            double W = Y - Math.Log(K);
            Complex y = Complex.Exp(i*u*W)*CFXT/(u*u + 0.25);

            // Return the real part of the integral
            return y.Real;
        }
        // The price in Equation 3.11 of Lewis (2001)
        public double LewisPrice311(double S,double K,double r,double q,double T,double theta,double kappa,double sigma,double rho,double v0,int trap,double[] x,double[] w)
        {
            double lambda = 0.0;
            Complex i = new Complex(0.0,1.0);

            // Compute the integral.
            // The integration variable is complex k = kr - (1/2)*i
            int Nx = x.Length;
            double[] Int = new double[Nx];
            Complex u = 0.0;
            for(int j=0;j<=Nx-1;j++)
            {
                u = x[j] - i/2.0;
                Int[j] = w[j] * LewisIntegrand311(u,kappa,theta,lambda,rho,sigma,T,S,K,r,q,v0,trap);
            }
            double Integral = Int.Sum();

            // Equation (3.11) in Lewis (2011)
            return S*Math.Exp(-q*T) - (1.0/Math.PI)*Math.Sqrt(K*S)*Math.Exp(-(r+q)*T/2.0)*Integral;
        }
    }
}

