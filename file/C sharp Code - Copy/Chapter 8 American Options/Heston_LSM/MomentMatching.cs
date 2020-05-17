using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_LSM
{
    class MMSimulation
    {
        // Simulation of stock price paths and variance paths using Euler or Milstein schemes
        public double[,] MMSim(HParam param,OpSet settings,int NT,int NS,double[,] Zv,double[,] Zs)
        {
            // Heston parameters
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            // option settings
            double r = settings.r;
            double q = settings.q;

            // Time increment
            double dt = settings.T/Convert.ToDouble(NT);

            // Initialize the variance and stock processes
            double[,] V = new double[NT,NS];
            double[,] S = new double[NT,NS];

            // Flags for negative variances
            int F = 0;

            // Starting values for the variance and stock processes
            for(int s=0;s<=NS-1;s++)
            {
                S[0,s] = settings.S;       // Spot price 
                V[0,s] = v0;       // Heston v0 initial variance 
            }

            // Generate the stock and volatility paths
            for(int s=0;s<=NS-1;s++)
            {
                for(int t=1;t<=NT-1;t++)
                {
                    // Matched moment lognormal approximation
                    double dW = Math.Sqrt(dt)*Zv[t,s];
                    double num = 0.5*sigma*sigma*V[t-1,s]*(1.0-Math.Exp(-2.0*kappa*dt)) / kappa;
                    double den = Math.Pow(Math.Exp(-kappa*dt)*V[t-1,s] + (1.0-Math.Exp(-kappa*dt))*theta,2);
                    double Gam = Math.Log(1 + num/den);
                    V[t,s] = (Math.Exp(-kappa*dt)*V[t-1,s] + (1.0-Math.Exp(-kappa*dt))*theta) * Math.Exp(-0.5*Gam*Gam + Gam*Zv[t,s]);

                    // Euler/Milstein discretization scheme for the log stock prices
                    S[t,s] = S[t-1,s]*Math.Exp((r-q-V[t-1,s]/2)*dt + Math.Sqrt(V[t-1,s]*dt)*Zs[t,s]);
                }
            }
            return S;
        }
    }
}


