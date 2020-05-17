using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Attari_Delta_and_Gamma
{
    class AttariGreeks
    {
        // Attari Integrand for the Greeks
        public double AttariProbGreeks(double u,double kappa,double theta,double lambda,double rho,double sigma,double T,
                          double K,double S0,double r,double q,double v0,int Trap,string Greek)
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

            // The L Function
            double L = Math.Log(Math.Exp(-r*T)*K/S0);

            // The coefficients
            double F = (f.Real + f.Imaginary/u);
            double H = (f.Imaginary - f.Real/u);
            
            // Return the integrand for the chosen Greek
            double A = 0.0;
            double y = 0.0;
            if (Greek == "Delta")
                A = (F*Math.Sin(L*u) - H*Math.Cos(L*u)) * u/S0/(1.0+u*u);
            else if(Greek == "Gamma")
            {
                y = (-F*(Math.Cos(L*u)*u + Math.Sin(L*u)) + H*(-Math.Sin(L*u)*u + Math.Cos(L*u)));
                A = y * u/S0/S0/(1.0+u*u);
            }
            return A;
        }
        // Analytic Attari Greeks
        public double AttariGreeksAnalytic(string PutCall,double kappa,double theta,double lambda,double rho,double sigma,double T,double K,double S0,double r,double q,double v0,int trap,string Greek,double[] x,double[] w)
        {
            int N = x.Length;
            double PI = Math.PI;
            double y = 0.0;
            double[] A = new double[N];
            for(int j=0;j<=N-1;j++)
            {
                if(Greek == "Delta")
                    A[j] = w[j] * AttariProbGreeks(x[j],kappa,theta,lambda,rho,sigma,T,K,S0,r,q,v0,trap,"Delta");
                else if(Greek == "Gamma")
                    A[j] = w[j] * AttariProbGreeks(x[j],kappa,theta,lambda,rho,sigma,T,K,S0,r,q,v0,trap,"Gamma");

                double Integral = A.Sum();

                if(Greek == "Gamma")
                    y = -Math.Exp(-r*T)*K/PI*Integral;
                else if(Greek == "Delta")
                    if(PutCall == "C")
                        y = 1.0 - Math.Exp(-r*T)*K/PI*Integral;
                    else
                        y = -Math.Exp(-r*T)*K/PI*Integral;
            }
            return y;
        }
    }
}
