using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Bivariate_Characteristic_Function
{
    class CharFun
    {
        // Heston Univariate characteristic function (f2)
        public Complex HestonCF(Complex phi,HParam param,OpSet settings)
        {
            Complex i     = new Complex(0.0,1.0);                   // Imaginary unit
            Complex S     = new Complex(settings.S,0.0);	                // Spot Price
            Complex K     = new Complex(settings.K,0.0);               // Strike Price
            Complex T     = new Complex(settings.T,0.0);               // Maturity in years
            Complex r     = new Complex(settings.r,0.0);               // Interest rate
            Complex q     = new Complex(settings.q,0.0);               // Dividend yield
            Complex rho   = new Complex(param.rho,0.0);         // Heston parameter: correlation
            Complex kappa = new Complex(param.kappa,0.0);         // Heston parameter: mean reversion speed
            Complex theta = new Complex(param.theta,0.0);         // Heston parameter: mean reversion speed
            Complex lambda = new Complex(param.lambda,0.0);         // Heston parameter: price of volatility risk
            Complex sigma = new Complex(param.sigma,0.0);         // Heston parameter: volatility of variance
            Complex v0    = new Complex(param.v0,0.0);         // Heston parameter: initial variance
            Complex x     = Complex.Log(S);
            Complex a     = kappa * theta;
            int Trap      = settings.trap;
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

            // The characteristic function.
            return Complex.Exp(C + D*v0 + i*phi*x);
        }
        public Complex HestonBivariateCF(Complex phi1,Complex phi2,HParam param,OpSet settings)
        {
            Complex i     = new Complex(0.0,1.0);               // Imaginary unit
            Complex S     = new Complex(settings.S,0.0);	    // Spot Price
            Complex K     = new Complex(settings.K,0.0);        // Strike Price
            Complex T     = new Complex(settings.T,0.0);        // Maturity in years
            Complex r     = new Complex(settings.r,0.0);        // Interest rate
            Complex q     = new Complex(settings.q,0.0);        // Dividend yield
            Complex rho   = new Complex(param.rho,0.0);         // Heston parameter: correlation
            Complex kappa = new Complex(param.kappa,0.0);       // Heston parameter: mean reversion speed
            Complex theta = new Complex(param.theta,0.0);       // Heston parameter: mean reversion speed
            Complex lambda = new Complex(param.lambda,0.0);     // Heston parameter: price of volatility risk
            Complex sigma = new Complex(param.sigma,0.0);       // Heston parameter: volatility of variance
            Complex v0    = new Complex(param.v0,0.0);          // Heston parameter: initial variance
            Complex x     = Complex.Log(S);
            Complex a     = kappa * theta;
            int trap      = settings.trap;
            Complex b,u,d,g,c,A=0.0,C=0.0,B,G = new Complex();

            // Log stock price
            Complex x0 = Complex.Log(S);

            u = -0.5;
            b = kappa + lambda;

            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi1 - b,2) - sigma*sigma*(2.0*u*i*phi1 - phi1*phi1));
            g = (b - rho*sigma*i*phi1 + d - sigma*sigma*i*phi2) 
              / (b - rho*sigma*i*phi1 - d - sigma*sigma*i*phi2);
            c = 1/g;
            B = i*phi1;

            if(trap==1)
            {
                // Little Trap formulation in Kahl (2008)
                G = (c*Complex.Exp(-d*T)-1.0)/(c-1.0);
                A = (r-q)*i*phi1*T + a/sigma/sigma*((b - rho*sigma*i*phi1 - d)*T - 2.0*Complex.Log(G));
                C = ((b - rho*sigma*i*phi1 - d) - (b - rho*sigma*i*phi1 + d)*c*Complex.Exp(-d*T)) / (sigma*sigma) / (1-c*Complex.Exp(-d*T));
            }
            else if(trap==0)
            {
                // Original Heston formulation.
                G = (1.0-g*Complex.Exp(d*T))/(1.0-g);
                A = (r-q)*i*phi1*T + a/sigma/sigma*((b - rho*sigma*i*phi1 + d)*T - 2.0*Complex.Log(G));
                C = ((b - rho*sigma*i*phi1 + d) - (b - rho*sigma*i*phi1 - d)*g*Complex.Exp(d*T)) / (sigma*sigma) / (1.0-g*Complex.Exp(d*T));
            }

            // The characteristic function.
            return Complex.Exp(A + B*x0 + C*v0);
        }
    }
}

