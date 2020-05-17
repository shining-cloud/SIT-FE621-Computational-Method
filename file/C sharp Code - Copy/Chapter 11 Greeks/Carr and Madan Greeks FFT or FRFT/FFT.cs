using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Carr_and_Madan_FFT_or_FRFT_Greeks
{
    class FastFT
    {
        // Fast Fourier Transform
        public double[,] HestonFFTGreek(int N,double uplimit,double alpha,string rule,HParam param,OpSet settings,string Greek)
        {
            HestonGreeks HG = new HestonGreeks();
            
            double s0 = Math.Log(settings.S);
            double pi = Math.PI;

            // Specify the increments
            double eta = uplimit/Convert.ToDouble(N);
            double lambdainc = 2*pi/Convert.ToDouble(N)/eta;

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

            double tau = settings.T;
            double r   = settings.r;
            double q   = settings.q;

            // Implement the FFT
            Complex i = new Complex(0.0,1.0);
            Complex[] psi = new Complex[N];
            Complex[] phi = new Complex[N];
            Complex[] x   = new Complex[N];
            Complex[] e   = new Complex[N];
            double[] CallFFT = new double[N];
            double[] sume = new double[N];
            for(int u=0;u<=N-1;u++)
            {
                for(int j=0;j<=N-1;j++)
                {
                    psi[j] = HG.HestonCFGreek(v[j]-(alpha+1.0)*i,param,settings,Greek);
                    phi[j] = Complex.Exp(-r*tau)*psi[j]/(alpha*alpha + alpha - v[j]*v[j] + i*v[j]*(2.0*alpha+1.0));
                    x[j]   = Complex.Exp(i*(b-s0)*v[j])*phi[j]*w[j];
                    e[j]   = Complex.Exp(-i*2.0*pi/Convert.ToDouble(N)*j*u)*x[j];
                    sume[u] += e[j].Real;
                }
                CallFFT[u] = eta*Math.Exp(-alpha*k[u])/pi * sume[u];
            }

            // Return the FFT call price and the strikes
            double[,] output = new double[N,2];
            for(int j=0;j<=N-1;j++)
            {
                output[j,0] = K[j];
                output[j,1] = CallFFT[j];
            }
            return output;
        }
    }
}


