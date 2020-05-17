using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Estimation_on_SP500_by_SVC
{
    class HestonPrice
    {
        // Heston characteristic function (f2) ==============================================================================
        public Complex HestonCF(Complex phi,HParam param,double S,double r,double q,double T,int trap)
        {
            Complex i = new Complex(0.0,1.0);                   // Imaginary unit
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0 = param.v0;
            double rho = param.rho;
            double lambda = param.lambda;
            double x = Math.Log(S);
            double a = kappa*theta;
            Complex b,u,d,g,c,D,G,C = new Complex();

            u = -0.5;
            b = kappa + lambda;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(trap==1)
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

        // Attari Characteristic Function ==================================================================================================
        public Complex AttariCF(double phi,HParam param,double Tau,double Spot,double Rate,double Div,int Trap)
        {
            Complex S     = new Complex(Spot,0.0);	    // Spot Price
            Complex T     = new Complex(Tau,0.0);       // Maturity in years
            Complex r     = new Complex(Rate,0.0);       // Interest rate
            Complex q     = new Complex(Div,0.0);       // Dividend yield
            Complex i     = new Complex(0.0,1.0);       // Imaginary unit
            Complex rho   = new Complex(param.rho,0.0);       // Heston parameter: correlation
            Complex kappa = new Complex(param.kappa,0.0);       // Heston parameter: mean reversion speed
            Complex theta = new Complex(param.theta,0.0);       // Heston parameter: mean reversion speed
            Complex lambda = new Complex(param.lambda,0.0);       // Heston parameter: price of volatility risk
            Complex sigma = new Complex(param.sigma,0.0);       // Heston parameter: volatility of variance
            Complex v0    = new Complex(param.v0,0.0);        // Heston parameter: initial variance
            Complex x     = Complex.Log(S);
            Complex a     = kappa * theta;
            Complex b,u,d,g,c,D,G,C = new Complex();

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

            // The Attari characteristic function.
            return Complex.Exp(C + D*v0 + i*phi*x) * Complex.Exp(-i*phi*(x+r*T));
        }

        // Heston Price by Gauss-Laguerre Integration =============================================================================================================
        public double HestonPriceGaussLaguerre(HParam param,double S,double K,double r,double q,double T,int trap,string PutCall,double[] X,double[] W)
        {
            Complex i = new Complex(0.0,1.0);
            Complex[] f1 = new Complex[32];
            Complex[] f2 = new Complex[32];
            double[] int1 = new double[32];
            double[] int2 = new double[32];

            // Numerical integration
            for(int j=0;j<=31;j++)
            {
                double phi = X[j];
                f1[j] = HestonCF(phi-i,param,S,r,q,T,trap) / (S*Math.Exp((r-q)*T));
                f2[j] = HestonCF(phi,param,S,r,q,T,trap);
                Complex FF1 = Complex.Exp(-i*phi*Complex.Log(K))*f1[j]/i/phi;
                Complex FF2 = Complex.Exp(-i*phi*Complex.Log(K))*f2[j]/i/phi;
                int1[j] = W[j] * FF1.Real;
                int2[j] = W[j] * FF2.Real;
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
    }
}


