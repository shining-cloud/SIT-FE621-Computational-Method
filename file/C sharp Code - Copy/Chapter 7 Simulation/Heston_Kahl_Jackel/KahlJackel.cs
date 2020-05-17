using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Kahl_Jackel
{
    class KahlJackel
    {
        // Price by simulation
        public double KahlJackelPrice(string scheme,string negvar,HParam param,OpSet settings,double alpha,int NT,int NS,string PutCall)
        {
            RandomNumbers RN = new RandomNumbers();
            double[] STe = KahlJackelSim(scheme,negvar,param,settings,alpha,NT,NS);
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
        public double[] KahlJackelSim(string scheme,string negvar,HParam param,OpSet settings,double alpha,int NT,int NS)
        {
            RandomNumbers RN = new RandomNumbers();

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
            double Zv,Zs;
            for(int s=0;s<=NS-1;s++)
            {
                for(int t=1;t<=NT-1;t++)
                {
                    // Generate two dependent N(0,1) variables with correlation rho
                    Zv = RN.RandomNorm();
                    Zs = rho*Zv + Math.Sqrt(1.0-rho*rho)*RN.RandomNorm();

                    if(scheme == "IJK")
                    {
                        // Implicit Milstein for the variance.
                        V[t,s] = (V[t-1,s] + kappa*theta*dt + sigma*Math.Sqrt(V[t-1,s]*dt)*Zv + sigma*sigma*dt*(Zv*Zv-1.0)/4.0) / (1+kappa*dt);

                        // Apply the full truncation or reflection scheme to the variance
                        if(V[t,s] <= 0.0)
                        {
                            F += 1;
                            if(negvar == "Reflection")          // Reflection: take -V
                                V[t,s] = Math.Abs(V[t,s]);
                            else if(negvar == "Truncation")
                                V[t,s] = Math.Max(0.0,V[t,s]);   // Truncation: take max(0,V)
                        }

                        // IJK discretization scheme for the log stock prices
                        S[t,s] = S[t-1,s]*Math.Exp((r-q-(V[t,s]+V[t-1,s])/4.0)*dt
            				   + rho*Math.Sqrt(V[t-1,s]*dt)*Zv
			                   + 0.5*(Math.Sqrt(V[t-1,s])+Math.Sqrt(V[t,s]))*(Zs-rho*Zv)*Math.Sqrt(dt)
            				   + rho*sigma*dt*(Zv*Zv-1.0)/2.0);
                    }
                    else if(scheme  == "PW")
                    {
                        // Pathwise Adapated Linearization Quadratic for the variance.
                        // Note: Bn = dZ/dt = Z*Math.Sqrt(dt)/dt = Z/Math.Sqrt(dt);
                        double theta2 = theta - sigma*sigma/4.0/kappa;
                        double Bn = Zv/Math.Sqrt(dt);
                        V[t,s] = V[t-1,s]
                               + (kappa*(theta2-V[t-1,s]) + sigma*Bn*Math.Sqrt(V[t-1,s]))*dt
                			   * (1.0 + (sigma*Bn-2.0*kappa*Math.Sqrt(V[t-1,s]))*dt/4.0/Math.Sqrt(V[t-1,s]));

                        // Apply the full truncation or reflection scheme to the variance
                        if(V[t,s] <= 0.0)
                        {
                            F += 1;
                            if(negvar == "Reflection")          // Reflection: take -V
                                V[t,s] = Math.Abs(V[t,s]);
                            else if(negvar == "Truncation")
                                V[t,s] = Math.Max(0.0,V[t,s]);   // Truncation: take max(0,V)
                        }

                        // Euler/Milstein discretization scheme for the log stock prices
                        S[t,s] = S[t-1,s]*Math.Exp((r-q-V[t-1,s]/2.0)*dt + Math.Sqrt(V[t-1,s]*dt)*Zs);
                    }
                    else if(scheme == "B")
                    {
                        // Balanced Implicit scheme for the variance
                        double absdW = Math.Sqrt(dt)*Math.Abs(Zv);
                        double C = kappa*dt + sigma/Math.Sqrt(V[t-1,s])*Math.Sqrt(dt)*Math.Abs(Zv);
                        V[t,s] = (V[t-1,s]*(1.0+C) + kappa*(theta-V[t-1,s])*dt + sigma*Math.Sqrt(V[t-1,s]*dt)*Zv) / (1.0+C);

                        // Euler/Milstein discretization scheme for the log stock prices
                        S[t,s] = S[t-1,s]*Math.Exp((r-q-V[t-1,s]/2.0)*dt + Math.Sqrt(V[t-1,s]*dt)*Zs);
                    }

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


