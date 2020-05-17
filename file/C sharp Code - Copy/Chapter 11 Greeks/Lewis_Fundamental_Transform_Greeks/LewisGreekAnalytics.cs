using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Fundamental_Transform_Greeks
{
    class LewisGreeks
    {
        // Lewis (2000) Integrand
        public double LewisIntegrandGreek(Complex k,double S,double r,double q,double X,double v0,double tau,double Theta,double Kappa,double sigma,double rho,string Greek)
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

            // The integrands
            Complex y = 0.0,dC,dD,I;
            if(Greek == "Price")
                y = Complex.Exp(-X*i*k-r*tau)/(k*k-i*k)*H;
            else if(Greek == "Delta")
                y = Complex.Exp(-X*i*k-r*tau)/(k*k-i*k)*H*(-i*k/S);
            else if(Greek == "Gamma")
                y = -Complex.Exp(-X*i*k-r*tau)*H/S/S;
            else if(Greek == "Vega1")
                y = Complex.Exp(-X*i*k-r*tau)/(k*k-i*k)*H*D*2.0*Math.Sqrt(v0);
            else if(Greek == "Rho")
                y = Complex.Exp(-X*i*k-r*tau)*(-i*k*tau-tau)/(k*k-i*k)*H;
            else if(Greek == "Theta")
            {
                dC = kappa*theta*((kappa+d)/2.0 + g*d*Complex.Exp(d*t)/(1.0-g*Complex.Exp(d*t))) * sigma*sigma/2.0;
                dD = (kappa+d)/2.0*Complex.Exp(d*t)*d*(g-1.0)/Complex.Pow(1-g*Complex.Exp(d*t),2.0) * sigma*sigma/2.0;
                y =  (-i*k*(r-q)-r + (dC+dD*v0))*H*Complex.Exp(-X*i*k-r*tau)/(k*k-i*k);
            }
            else if(Greek == "Volga")
            {
                I = Complex.Exp(-X*i*k-r*tau)/(k*k-i*k)*H*D*D*Math.Sqrt(v0)
                  + Complex.Exp(-X*i*k-r*tau)/(k*k-i*k)*H*D/2.0/Math.Sqrt(v0);
                y = I*4.0*Math.Sqrt(v0);
            }
            else if(Greek == "Vanna")
                y = Complex.Exp(-X*i*k-r*tau)/(k*k-i*k)*H*D*(-i*k/S)*2.0*Math.Sqrt(v0);

            // Return the real part of the integrand.
            return y.Real;
        }

        // The Lewis Call price by Gauss-Laguerre Integration
        public double HestonLewisGreekPrice(double S,double K,double r,double q,double v0,double tau,double ki,
                                           double theta,double kappa,double sigma,double rho,int form,
                                           double[] x,double[] w,string Greek)
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
            int Nx = x.Length;
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            Complex u;
            double[] int1 = new Double[Nx];
            for(int k=0;k<=Nx-1;k++)
            {
                u = x[k] + i*ki;
                int1[k] = w[k] * LewisIntegrandGreek(u,S,r,q,X,v0,tau,theta,kappa,sigma,rho,Greek);
            }
            double Integral = int1.Sum();

            // Compute the Greeks
            double pi = Math.PI;
            switch(Greek)
            {
                case "Price":
                    if(form==2)
                        // Price based on C2(K) Equation (2.10) page 41 of Lewis (2000)
                        return S*Math.Exp(-q*tau) - (K/pi)*Integral;
                    else
                        // Price based on C1(K) Equation (2.8) page 40 of Lewis (2000)
                        return -(K/pi)*Integral;
                    break;
                case "Delta":
                    if(form==2)
                        // Delta based on C2(K) Equation (2.10) page 41 of Lewis (2000)
                        return Math.Exp(-q*tau) - (K/pi)*Integral;
                    else
                        // Delta based on C1(K) Equation (2.8) page 40 of Lewis (2000)
                        return -(K/pi)*Integral;
                    break;
                case "Theta":
                    if(form==2)
                        return q*S*Math.Exp(-q*tau) + (K/pi)*Integral;
                    else
                        return (K/pi)*Integral;
                default:
                    // Remaining Greeks based on C1(K) or C2(K) have the same integrand
                    return -(K/pi)*Integral;
                    break;
            }
        }
    }
}




