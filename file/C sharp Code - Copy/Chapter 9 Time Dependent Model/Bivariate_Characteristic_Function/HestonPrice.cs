using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Bivariate_Characteristic_Function
{
    class HestonPrice
    {   
        // Heston Price by Gauss-Laguerre Integration
        public double HestonPriceGaussLaguerre(HParam param, OpSet settings,double[] X, double[] W, int CF)
        {
            Complex i = new Complex(0.0,1.0);
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;
            double S = settings.S;
            double K = settings.K;
            double r = settings.r;
            double q = settings.q;
            double T = settings.T;
            int trap = settings.trap;

            int N = X.Length;
            Complex f1,f2,I1,I2;
            double[] int1 = new Double[N];
            double[] int2 = new Double[N];
            double phi;

            // Numerical integration
            CharFun HCF = new CharFun();
            for(int k=0;k<=31;k++)
            {
                f2 = 0;
                f1 = 0;
                phi = X[k];
                if(CF==1)
                {
                    f2 = HCF.HestonCF(phi,param,settings);
                    f1 = HCF.HestonCF(phi-i,param,settings) / (S*Math.Exp((r-q)*T));
                }
                else if(CF==2)
                {
                    f2 = HCF.HestonBivariateCF(phi,  0.0,param,settings);
                    f1 = HCF.HestonBivariateCF(phi-i,0.0,param,settings) / (S*Math.Exp((r-q)*T));
                }
                I2 = Complex.Exp(-i*phi*Complex.Log(K))*f2/i/phi;
                I1 = Complex.Exp(-i*phi*Complex.Log(K))*f1/i/phi;
                int2[k] = W[k] * I2.Real;
                int1[k] = W[k] * I1.Real;
            }

            // Define P1 and P2
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*int1.Sum();
            double P2 = 0.5 + 1.0/pi*int2.Sum();

            // The call price
            string PutCall = settings.PutCall;
            double HestonC = S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);

            // Output the option price
            if (PutCall == "C")
	            return HestonC;
            else
	            return HestonP;
 
        }
    }
}

