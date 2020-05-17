using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Moment_Matching
{
    class MMSimulation
    {
        // Price by simulation
        public double MMPrice(HParam param,OpSet settings,int NT,int NS)
        {
            RandomNumbers RN = new RandomNumbers();
            double[] STe = MMSim(param,settings,NT,NS);
            double[] Price = new double[NS];
            for(int s=0;s<=NS-1;s++)
            {
                if(settings.PutCall == "C")
                    Price[s] = Math.Max(STe[s] - settings.K,0.0);
                else if(settings.PutCall == "P")
                    Price[s] = Math.Max(settings.K - STe[s],0.0);
            }
            return Math.Exp(-settings.r*settings.T) * RN.Mean(Price);
        }

        // Simulation of stock price paths and variance paths using Euler or Milstein schemes
        public double[] MMSim(HParam param,OpSet settings,int NT,int NS)
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
            RandomNumbers RN = new RandomNumbers();
            double Zv,Zs;
            for(int s=0;s<=NS-1;s++)
            {
                for(int t=1;t<=NT-1;t++)
                {
                    // Generate two dependent N(0,1) variables with correlation rho
                    Zv = RN.RandomNorm();
                    Zs = rho*Zv + Math.Sqrt(1.0-rho*rho)*RN.RandomNorm();

                    // Matched moment lognormal approximation
                    double dW = Math.Sqrt(dt)*Zv;
                    double num = 0.5*sigma*sigma*V[t-1,s]*(1.0-Math.Exp(-2.0*kappa*dt)) / kappa;
                    double den = Math.Pow(Math.Exp(-kappa*dt)*V[t-1,s] + (1.0-Math.Exp(-kappa*dt))*theta,2);
                    double Gam = Math.Log(1.0 + num/den);
                    V[t,s] = (Math.Exp(-kappa*dt)*V[t-1,s] + (1.0-Math.Exp(-kappa*dt))*theta) * Math.Exp(-0.5*Gam*Gam + Gam*Zv);

                    // Euler/Milstein discretization scheme for the log stock prices
                    S[t,s] = S[t-1,s]*Math.Exp((r-q-V[t-1,s]/2.0)*dt + Math.Sqrt(V[t-1,s]*dt)*Zs);
                }
            }
            // Return the vector or terminal stock prices
            double[] output = new double[NS];
            for(int s=0;s<=NS-1;s++)
                output[s] = S[NT-1,s];

            return output;
        }
    }
}


