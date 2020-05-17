using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Beliaeva_Nawalkha_Tree
{
    partial class BNTree
    {
        // Heston Integrand
        static float HestonProb(float phi,HParam param,float S,float K,float r,float q,float T,int Pnum,int Trap)
        {
            Complex i  = new Complex(0.0,1.0);                   // Imaginary unit
            float kappa = param.kappa;
            float theta = param.theta;
            float sigma = param.sigma;
            float v0 = param.v0;
            float rho = param.rho;
            float lambda = 0.0f;
            float x = Convert.ToSingle(Math.Log(S));
            float a = kappa*theta;
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

            // Return the real part of the integrand.
            return Convert.ToSingle(integrand.Real);
        }

        // Heston Price by Gauss-Laguerre Integration =================================================================
        static float HestonPriceGaussLaguerre(HParam param,float S,float K,float r,float q,float T,int trap,string PutCall,float[] x,float[] w)
        {
            float[] int1 = new float[32];
            float[] int2 = new float[32];

            // Numerical integration
            for(int j=0;j<=31;j++)
            {
                int1[j] = w[j] * HestonProb(x[j],param,S,K,r,q,T,1,trap);
                int2[j] = w[j] * HestonProb(x[j],param,S,K,r,q,T,2,trap);
            }

            // Define P1 and P2
            float pi = Convert.ToSingle(Math.PI);
            float P1 = 0.5f + 1.0f/pi*int1.Sum();
            float P2 = 0.5f + 1.0f/pi*int2.Sum();

            // The call price
            float HestonC = Convert.ToSingle(S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2);

            // The put price by put-call parity
            float HestonP = Convert.ToSingle(HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T));

            // Output the option price
            if(PutCall == "C")
                return HestonC;
            else
                return HestonP;
        }
    }
}

