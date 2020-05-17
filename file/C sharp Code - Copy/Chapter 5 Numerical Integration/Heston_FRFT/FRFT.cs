using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_FRFT
{
    class HFRFT
    {
        // Call price by Fractional Fast Fourier Transform =====================================================================
        public double[,] HestonFRFT(int N,double uplimit,double alpha,string rule,double lambdainc,double eta,HParam param,OpSet settings)
        {
            double s0 = Math.Log(settings.S);
            double pi = Math.PI;

            // Initialize and specify the weights
            double[] w = new double[N];
            if(rule == "Trapezoidal")
            {
                w[0] = 0.5;
                w[N-1] = 0.5;
                for(int j=1;j<=N-2;j++)
                    w[j] = 1;
            }
            else if(rule == "Simpsons")
            {
                w[0] = 1.0/3.0;
                w[N-1] = 1.0/3.0;
                for(int j=1;j<=N-1;j++)
                    w[j] = (1.0/3.0) * (3 + Math.Pow(-1,j+1));
            }

            // Specify the b parameter
            double b = Convert.ToDouble(N)*lambdainc/2.0;

            // Create the grid for the integration
            double[] v = new double[N];
            for(int j=0;j<=N-1;j++)
                v[j] = eta * j;

            // Create the grid for the log-strikes and strikes
            double[] k = new double[N];
            double[] K = new double[N];
            for(int j=0;j<=N-1;j++)
            {
                k[j] = -b + lambdainc*Convert.ToDouble(j) + s0;
                K[j] = Math.Exp(k[j]);
            }

            // Beta parameter for the FRFT
            double beta = lambdainc*eta/2.0/pi;

            // Option price settings
            double tau = settings.T;
            double r   = settings.r;
            double q   = settings.q;

            // Implement the FFT
            Complex i = new Complex(0.0,1.0);
            Complex[] psi = new Complex[N];
            Complex[] phi = new Complex[N];
            Complex[] x   = new Complex[N];
            Complex[] y   = new Complex[N];
            double[] CallFRFT = new double[N];
            double[] sume = new double[N];

            // Build the input vectors for the FRFT
            HestonPrice HP = new HestonPrice();
            for(int j=0;j<=N-1;j++)
            {
                psi[j] = HP.HestonCF(v[j]-(alpha+1.0)*i,param,settings);
                phi[j] = Complex.Exp(-r*tau)*psi[j] / (alpha*alpha + alpha - v[j]*v[j] + i*v[j]*(2.0*alpha+1.0));
                x[j]   = Complex.Exp(i*(b-s0)*v[j]) * phi[j] * w[j];
            }
            y = FRFT(x,beta);
            Complex Call = new Complex(0.0,0.0);
            for(int u=0;u<=N-1;u++)
            {
                Call = eta*Complex.Exp(-alpha*k[u])*y[u]/pi;
                CallFRFT[u] = Call.Real;
            }

            // Return the FRFT call price and the strikes
            double[,] output = new double[N,2];
            for(int j=0;j<=N-1;j++)
            {
                output[j,0] = K[j];
                output[j,1] = CallFRFT[j];
            }
            return output;
        }
        // Fractional Fast Fourier Transform ======================================================================
        public Complex[] FRFT(Complex[] x,double beta)
        {
            int N = x.Length;
            Complex i = new Complex(0.0,1.0);
            double pi = Math.PI;

            // Construct the y and z vectors
            Complex[] y = new Complex[2*N];
            Complex[] z = new Complex[2*N];
            for(int j=0;j<=N-1;j++)
            {
                double J = Convert.ToDouble(j);
                y[j] = Complex.Exp(-i*pi*J*J*beta) * x[j];
                z[j] = Complex.Exp(i*pi*J*J*beta);
            }
            for(int j=N;j<=2*N-1;j++)
            {
                y[j] = 0.0;
                double M = Convert.ToDouble(2*N-j);
                z[j] = Complex.Exp(i*pi*M*M*beta);
            }

            // FFT on y and z
            Complex[] Dy = new Complex[2*N];
            Complex[] Dz = new Complex[2*N];
            Dy = FFT(y);
            Dz = FFT(z);

            // h vectors
            Complex[] h  = new Complex[2*N];
            for(int j=0;j<=2*N-1;j++)
                h[j] = Dy[j]*Dz[j];
            Complex[] ih = new Complex[2*N];
            ih = IFFT(h);

            // e vector
            Complex[] e = new Complex[N];
            for(int j=0;j<=N-1;j++)
            {
                double J = Convert.ToDouble(j);
                e[j] = Complex.Exp(-i*pi*J*J*beta);
            }

            // The FRFT vector of size N
            Complex[] xhat = new Complex[N];
            for(int j=0;j<=N-1;j++)
                xhat[j] = e[j] * ih[j];

            // Return the FRFT vector 
            return xhat;
        }

        // Fast Fourier Transform ============================================================================
        public Complex[] FFT(Complex[] x)
        {
            int N = x.Length;
            double pi = Math.PI;
            Complex i = new Complex(0.0,1.0);
            Complex coeff = new Complex(0.0,0.0);
            Complex[] xhat = new Complex[N];

            for(int k=0;k<=N-1;k++)
            {
                coeff = 0.0;
                for(int j=0;j<=N-1;j++)
                {
                    double K = Convert.ToDouble(k);
                    double J = Convert.ToDouble(j);
                    double M = Convert.ToDouble(N);
                    coeff += Complex.Exp(-i*2*pi*J*K/M) * x[j];
                }
                xhat[k] = coeff;
            }
            return xhat;
        }
        // Inverse Fast Fourier Transform ====================================================================
        public Complex[] IFFT(Complex[] xhat)
        {
            int N = xhat.Length;
            double pi = Math.PI;
            Complex i = new Complex(0.0,1.0);
            Complex coeff = new Complex(0.0,0.0);
            Complex[] x = new Complex[N];

            for(int k=0;k<=N-1;k++)
            {
                coeff = 0.0;
                for(int j=0;j<=N-1;j++)
                {
                    double K = Convert.ToDouble(k);
                    double J = Convert.ToDouble(j);
                    double M = Convert.ToDouble(N);
                    coeff += Complex.Exp(i*2*pi*J*K/M) * xhat[j];
                }
                x[k] = coeff / Convert.ToDouble(N);
            }
            return x;
        }

    }
}


