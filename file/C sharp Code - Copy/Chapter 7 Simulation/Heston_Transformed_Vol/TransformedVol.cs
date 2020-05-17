using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Transformed_Volatility
{
    class TVSimulation
    {
        // Price by simulation
        public double TransVolPrice(string scheme,HParam param,OpSet settings,int NT,int NS)
        {
            RandomNumbers RN = new RandomNumbers();
            double[] STe = TVSim(scheme,param,settings,NT,NS);
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
        public double[] TVSim(string scheme,HParam param,OpSet settings,int NT,int NS)
        {
            // Heston parameters
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double sv0   = Math.Sqrt(v0);
            double rho   = param.rho;
            double lambda = param.lambda;

            // option settings
            double r = settings.r;
            double q = settings.q;

            // Time increment
            double dt = settings.T/Convert.ToDouble(NT);

            // Initialize the variance, volatility, stock processes, and log stock process
            double[,] w = new double[NT,NS];
            double[,] v = new double[NT,NS];
            double[,] S = new double[NT,NS];
            double[,] X = new double[NT,NS];

            // Starting values for the variance and stock processes
            for(int s=0;s<=NS-1;s++)
            {
                S[0,s] = settings.S;            // Spot price 
                X[0,s] = Math.Log(settings.S);  // Spot log price
                v[0,s] = v0;                    // Heston v0 initial variance
                w[0,s] = sv0;                   // Heston sqrt(v0) initial volatility
            }

            // Generate the stock and volatility paths
            RandomNumbers RN = new RandomNumbers();
            double Zv,Zx;
            double m1,m2,beta,thetav;
            for(int s=0;s<=NS-1;s++)
            {
                for(int t=1;t<=NT-1;t++)
                {
                    // Generate two dependent N(0,1) variables with correlation rho
                    Zv = RN.RandomNorm();
                    Zx = rho*Zv + Math.Sqrt(1.0-rho*rho)*RN.RandomNorm();

                    if(scheme == "Euler")
                        // Euler scheme for the volatility
                        w[t,s] = w[t-1,s] + 0.5*kappa*((theta-sigma*sigma/4.0/kappa)/w[t-1,s] - w[t-1,s])*dt + 0.5*sigma*Math.Sqrt(dt)*Zv;

                    else if(scheme =="TV")
                    {
                        // Zhu (2010) scheme
                        m1 = theta + (v[t-1,s] - theta)*Math.Exp(-kappa*dt);
                        m2 = sigma*sigma/4.0/kappa*(1.0-Math.Exp(-kappa*dt));
                        beta = Math.Sqrt(Math.Max(0.0,m1-m2));
                        thetav = (beta - w[t-1,s]*Math.Exp(-kappa*dt/2.0))/(1.0-Math.Exp(-kappa*dt/2.0));
                        w[t,s] = w[t-1,s] + 0.5*kappa*(thetav - w[t-1,s])*dt + 0.5*sigma*Math.Sqrt(dt)*Zv;
                    }
                    v[t,s] = w[t,s]*w[t,s];

                    // Discretize the log stock price
                    X[t,s] = X[t-1,s] + (settings.r-settings.q-w[t-1,s]*w[t-1,s]/2.0)*dt + w[t-1,s]*Math.Sqrt(dt)*Zx;
                    S[t,s] = Math.Exp(X[t,s]);
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

