using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Double_Heston_Simulation
{
    class HestonPriceDH
    {
        // Heston Integrand
        public Complex HestonDoubleProb(Complex phi,DHParam param,OpSet settings)
        {
            Complex i      = new Complex(0.0,1.0);                   // Imaginary unit
            double S = settings.S;
            double K = settings.K;
            double T = settings.T;
            double r = settings.r;
            double q = settings.q;
            double kappa1 = param.kappa1;
            double theta1 = param.theta1;
            double sigma1 = param.sigma1;
            double v01 = param.v01;
            double rho1 = param.rho1;
            double kappa2 = param.kappa2;
            double theta2 = param.theta2;
            double sigma2 = param.sigma2;
            double v02 = param.v02;
            double rho2 = param.rho2;
            double x0 = Math.Log(S);
            int Trap       = settings.trap;
            Complex d1,d2,G1,G2,B1,B2,X1,X2,A,g1,g2 = new Complex();

            if(Trap==1)
            {
                // Little Trap formulation of the characteristic function
                d1 = Complex.Sqrt(Complex.Pow(kappa1-rho1*sigma1*i*phi,2) + sigma1*sigma1*phi*(phi+i));
                d2 = Complex.Sqrt(Complex.Pow(kappa2-rho2*sigma2*i*phi,2) + sigma2*sigma2*phi*(phi+i));
                G1 = (kappa1-rho1*sigma1*phi*i-d1) / (kappa1-rho1*sigma1*phi*i+d1);
                G2 = (kappa2-rho2*sigma2*phi*i-d2) / (kappa2-rho2*sigma2*phi*i+d2);
                B1 = (kappa1-rho1*sigma1*phi*i-d1)*(1.0-Complex.Exp(-d1*T)) / (sigma1*sigma1) / (1.0-G1*Complex.Exp(-d1*T));
                B2 = (kappa2-rho2*sigma2*phi*i-d2)*(1.0-Complex.Exp(-d2*T)) / (sigma2*sigma2) / (1.0-G2*Complex.Exp(-d2*T));
                X1 = (1.0-G1*Complex.Exp(-d1*T))/(1.0-G1);
                X2 = (1.0-G2*Complex.Exp(-d2*T))/(1.0-G2);
                A  = (r-q)*phi*i*T 
		            + kappa1*theta1/sigma1/sigma1*((kappa1-rho1*sigma1*phi*i-d1)*T - 2.0*Complex.Log(X1))
		            + kappa2*theta2/sigma2/sigma2*((kappa2-rho2*sigma2*phi*i-d2)*T - 2.0*Complex.Log(X2));
            }
            else
            {
                // Original Heston formulation of the characteristic function
                d1 = Complex.Sqrt(Complex.Pow(kappa1-rho1*sigma1*phi*i,2) + sigma1*sigma1*(phi*i+phi*phi));
                d2 = Complex.Sqrt(Complex.Pow(kappa2-rho2*sigma2*phi*i,2) + sigma2*sigma2*(phi*i+phi*phi));
                g1 = (kappa1-rho1*sigma1*phi*i+d1)/(kappa1-rho1*sigma1*phi*i-d1);
                g2 = (kappa2-rho2*sigma2*phi*i+d2)/(kappa2-rho2*sigma2*phi*i-d2);
                B1 = (kappa1-rho1*sigma1*phi*i+d1)*(1.0-Complex.Exp(d1*T)) / (sigma1*sigma1) / (1.0-g1*Complex.Exp(d1*T));
                B2 = (kappa2-rho2*sigma2*phi*i+d2)*(1.0-Complex.Exp(d2*T)) / (sigma2*sigma2) / (1.0-g2*Complex.Exp(d2*T));
                X1 = (1.0-g1*Complex.Exp(d1*T))/(1.0-g1);
                X2 = (1.0-g2*Complex.Exp(d2*T))/(1.0-g2);
                A  = (r-q)*phi*i*T
		           + kappa1*theta1/sigma1/sigma1*((kappa1-rho1*sigma1*phi*i+d1)*T - 2.0*Complex.Log(X1))
		           + kappa2*theta2/sigma2/sigma2*((kappa2-rho2*sigma2*phi*i+d2)*T - 2.0*Complex.Log(X2));
            }

            // The characteristic function.
            return Complex.Exp(A + B1*v01 + B2*v02 + i*phi*x0);
        }

        // Heston Price by Trapezoidal rule
        public double DoubleHestonPriceNewtonCoates(DHParam param,OpSet settings,double a,double b,int N)
        {
            double S = settings.S;
            double K = settings.K;
            double r = settings.r;
            double q = settings.q;
            double T = settings.T;
            string PutCall = settings.PutCall;
            Complex i = new Complex(0.0,0.0);
            double h = (b-a)/Convert.ToDouble(N-1);
            double[] phi = new double[N];
            double[] w = new double[N];
            w[0] = 0.5*h;
            w[N-1] = 0.5*h;
            for(int k=0;k<=N-1;k++)
            {
                phi[k] = a + k*h;
                if(k>=1 & k<=N-2)
                    w[k] = h;
            }
            Complex u = new Complex(0.0,0.0);
            Complex[] f2 = new Complex[N];
            Complex[] f1 = new Complex[N];
            double[] Int1 = new double[N];
            double[] Int2 = new double[N];
            Complex I1,I2;
            for(int k=0;k<=N-1;k++)
            {
                u = phi[k];
                f2[k] = HestonDoubleProb(u,param,settings);
                f1[k] = HestonDoubleProb(u-i,param,settings);
                I2 = Complex.Exp(-i*u*Complex.Log(K))*f2[k]/i/u;
                I1 = Complex.Exp(-i*u*Complex.Log(K))*f1[k]/i/u/S/Complex.Exp((r-q)*T);
                Int2[k] = w[k] * I2.Real;
                Int1[k] = w[k] * I1.Real;
            }

            // The ITM probabilities
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*Int1.Sum();
            double P2 = 0.5 + 1.0/pi*Int2.Sum();

            // The Double Heston call price
            double HestonC = S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*T) + K*Math.Exp(-r*T);

            // Output the option price
            if(PutCall == "C")
                return HestonC;
            else
                return HestonP;
        }
        public double DoubleHestonPriceGaussLaguerre(DHParam param,OpSet settings,double[] X,double[] W)
        {
            int N = X.Length;
            double S = settings.S;
            double K = settings.K;
            double r = settings.r;
            double q = settings.q;
            double T = settings.T;
            string PutCall = settings.PutCall;
            Complex i = new Complex(0.0,1.0);
            Complex u = new Complex(0.0,0.0);
            Complex[] f2 = new Complex[N];
            Complex[] f1 = new Complex[N];
            double[] Int1 = new double[N];
            double[] Int2 = new double[N];
            Complex I1 = new Complex(0.0,0.0);
            Complex I2 = new Complex(0.0,0.0);
            for(int k=0;k<=N-1;k++)
            {
                u = X[k];
                f2[k] = HestonDoubleProb(u,param,settings);
                f1[k] = HestonDoubleProb(u-i,param,settings);
                I2 = Complex.Exp(-i*u*Complex.Log(K))*f2[k]/i/u;
                I1 = Complex.Exp(-i*u*Complex.Log(K))*f1[k]/i/u/S/Complex.Exp((r-q)*T);
                Int2[k] = W[k] * I2.Real;
                Int1[k] = W[k] * I1.Real;
            }

            // The ITM probabilities
            double pi = Math.PI;
            double P1 = 0.5 + 1.0/pi*Int1.Sum();
            double P2 = 0.5 + 1.0/pi*Int2.Sum();

            // The Double Heston call price
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

