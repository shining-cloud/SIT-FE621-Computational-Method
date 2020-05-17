using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Fundamental_Transform
{
    class LewisPrice
    {
        // Lewis (2000) Integrand
        public double LewisIntegrand(Complex k,double X,double v0,double tau,double Theta,double Kappa,double sigma,double rho)
        {
            // Define the updated "tilde" parameters
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            Complex kappa = 2.0*(Kappa + i*k*rho*sigma)/(sigma*sigma);
            Complex theta = kappa*Theta/(kappa + i*k*rho*sigma);
            Double t = tau*sigma*sigma/2.0;
            Complex c = (k*k-i*k)/(sigma*sigma);

            // Compute the required quantities
            Complex d = Complex.Sqrt(kappa*kappa + 4.0*c);
            Complex alpha = (-kappa + d)/2.0;
            Complex beta  = (-kappa - d)/2.0;
            Complex g = beta/alpha;

            // The functions B(t) and A(t)
            Complex D = (kappa+d)*(1.0-Complex.Exp(d*t)) / (1.0 - g*Complex.Exp(d*t))/2.0;
            Complex C = kappa*theta*((kappa+d)*t/2.0 - Complex.Log((1-g*Complex.Exp(d*t))/(1.0-g)));

            // The fundamental transform, H(k,S,t)
            Complex H = Complex.Exp(C + D*v0);

            // The integrand.
            Complex y = Complex.Exp(-X*i*k)/(k*k - i*k)*H;

            // Return the real part of the integrand.
            return y.Real;
        }

        // The Lewis Call price by Gauss-Laguerre Integration
        public double HestonLewisCallPrice(double S,double K,double r,double q,double v0,double tau,double ki,
                                           double theta,double kappa,double sigma,double rho,int form,
                                           double[] x,double[] w)
        {
            // INPUTS
            //   S = Spot price
            //   K = Strike price
            //   r = Risk free rate
            //   q = Dividend Yield
            //   v0 = Heston parameter, initial variance
            //   tau = Parameter
            //   ki = Im(k) strip cutoff
            //   theta = Heston parameter, mean reversion level
            //   kappa = Heston parameter, mean reversion speed
            //   sigma = Heston parameter, volatility of variance
            //   form  = 2 produces C2(S,K,t) requires 0 < Im[k] < 1
            //         = 1 produces C1(S,K,t) requires 1 < Im[k] < B
            //   x = Gauss Laguerre rule: Abscissas
            //   w = Gauss Laguerre rule: Weights

            // Lewis Parameters
            double kmax = Math.Floor(Math.Max(1000,10/Math.Sqrt(v0*tau)));
            double X = Math.Log(S/K) + (r-q)*tau;

            // Compute the integral using 32-point Gauss Laguerre
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            Complex u;
            double[] int1 = new Double[32];
            for(int k=0;k<=31;k++)
            {
                u = x[k] + i*ki;
                int1[k] = w[k] * LewisIntegrand(u,X,v0,tau,theta,kappa,sigma,rho);
            }
            double Integral = int1.Sum();

            // Compute the call price
            double pi = Math.PI;
            if(form==2)
                // C2(S,K,t) Equation (2.10) page 41 of Lewis (2000)
                return S*Math.Exp(-q*tau) - (1.0/pi)*K*Math.Exp(-r*tau)*Integral;
            else
                // C1(S,K,t) Equation (2.8) page 40 of Lewis (2000)
                return -(1.0/pi)*K*Math.Exp(-r*tau)*Integral;
        }
    }
}


