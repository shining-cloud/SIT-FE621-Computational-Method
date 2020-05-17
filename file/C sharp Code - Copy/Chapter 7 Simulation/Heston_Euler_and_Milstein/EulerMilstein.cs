using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Euler_and_Milstein
{
    class Simulation
    {
        // Price by simulation
        public double EulerMilsteinPrice(string scheme,string negvar,HParam param,OpSet settings,double alpha,int NT,int NS,string PutCall)
        {
            RandomNumbers RA = new RandomNumbers();
            double[] STe = EulerMilsteinSim(scheme,negvar,param,settings,alpha,NT,NS);
            double[] Price = new double[NS];
            for(int s=0;s<=NS-1;s++)
            {
                if(settings.PutCall == "C")
                    Price[s] = Math.Max(STe[s] - settings.K,0.0);
                else if(settings.PutCall == "P")
                    Price[s] = Math.Max(settings.K - STe[s],0.0);
            }
            return Math.Exp(-settings.r*settings.T) * RA.Mean(Price);
        }

        // Simulation of stock price paths and variance paths using Euler or Milstein schemes
        public double[] EulerMilsteinSim(string scheme,string negvar,HParam param,OpSet settings,double alpha,int NT,int NS)
        {
            RandomNumbers RN = new RandomNumbers();

            // Heston parameters
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

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
            double Zv,Zs;
            for(int s=0;s<=NS-1;s++)
            {
                for(int t=1;t<=NT-1;t++)
                {
                    // Generate two dependent N(0,1) variables with correlation rho
                    Zv = RN.RandomNorm();
                    Zs = rho*Zv + Math.Sqrt(1.0-rho*rho)*RN.RandomNorm();

                    if(scheme == "Euler")
                        // Euler discretization for the variance
                        V[t,s] = V[t-1,s] + kappa*(theta-V[t-1,s])*dt + sigma*Math.Sqrt(V[t-1,s]*dt)*Zv;
                    else if(scheme == "Milstein")
                        // Milstein discretization for the variance.
                        V[t,s] = V[t-1,s] + kappa*(theta-V[t-1,s])*dt + sigma*Math.Sqrt(V[t-1,s]*dt)*Zv + 0.25*sigma*sigma*dt*(Zv*Zv-1.0);
                    else if(scheme == "IM")
                        // Implicit Milstein for the variance.
                        V[t,s] = (V[t-1,s] + kappa*theta*dt + sigma*Math.Sqrt(V[t-1,s]*dt)*Zv + sigma*sigma*dt*(Zv*Zv-1.0)/4.0) / (1.0+kappa*dt);
                    else if(scheme == "WM")
                        // Weighted Explicit-Implicit Milstein Scheme
                        V[t,s] = (V[t-1,s] + kappa*(theta-alpha*V[t-1,s])*dt + sigma*Math.Sqrt(V[t-1,s]*dt)*Zv + sigma*sigma*dt*(Zv*Zv-1.0)/4.0) / (1.0+(1.0-alpha)*kappa*dt);

                    // Apply the full truncation or reflection scheme to the variance
                    if(V[t,s] <= 0)
                    {
                        F += 1;
                        if(negvar == "Reflection")          // Reflection: take -V
                            V[t,s] = Math.Abs(V[t,s]);
                        else if(negvar == "Truncation")
                            V[t,s] = Math.Max(0,V[t,s]);   // Truncation: take max(0,V)
                    }

                    // Discretize the log stock price
                    S[t,s] = S[t-1,s] * Math.Exp((settings.r-settings.q-V[t-1,s]/2.0)*dt + Math.Sqrt(V[t-1,s]*dt)*Zs);
                }
            }

            // Return the vector or terminal stock prices
            double[] output = new double[NS];
            for (int s=0; s<=NS-1; s++)
                output[s] = S[NT-1,s];

            return output;
        }
    }
}

