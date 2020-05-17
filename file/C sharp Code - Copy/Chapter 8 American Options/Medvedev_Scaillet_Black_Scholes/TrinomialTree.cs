using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Black_Scholes
{
    class TrinomialPrice
    {
        public double TrinomialTree(double Spot,double K,double r,double q,double T,double v,int N,string EuroAmer,string PutCall)
        {
            // Trinomial tree parameters and probabilities.
            double dt = T/Convert.ToDouble(N);
            double u = Math.Exp(v*Math.Sqrt(2.0*dt));
            double d = 1/u;
            double pu = Math.Pow(Math.Exp((r-q)*dt/2.0) - Math.Exp(-v*Math.Sqrt(dt/2.0)),2.0)/Math.Pow(Math.Exp(v*Math.Sqrt(dt/2.0)) - Math.Exp(-v*Math.Sqrt(dt/2)),2.0);
            double pd = Math.Pow(Math.Exp(v*Math.Sqrt(dt/2.0)) - Math.Exp((r-q)*dt/2),2.0)/Math.Pow(Math.Exp(v*Math.Sqrt(dt/2.0)) - Math.Exp(-v*Math.Sqrt(dt/2.0)),2.0);
            double pm = 1 - pu - pd;

            // Initialize the stock prices
            double[,] S = new double[2*N+1,N+1];
            S[N,0] = Spot;

            // Build the trinomial tree
            for(int j=0;j<=N;j++)
                for(int i=0;i<=2*j;i++)
                    S[i,j] = Spot*Math.Pow(u,j-i);

            // Compute terminal payoffs
            double[,] V = new double[2*N+1,N+1];
            for(int i=0;i<=2*N;i++)
            {
                if(PutCall == "C")
                    V[i,N] = Math.Max(S[i,N] - K,0.0);
                else
                    V[i,N] = Math.Max(K - S[i,N],0.0);
            }

            // Backward recursion through the tree
            for(int j=N-1;j>=0;j--)
                for(int i=0;i<=2*j;i++)
                {
                    if(EuroAmer == "E")
                        V[i,j] = Math.Exp(-r*dt)*(pu*V[i,j+1] + pm*V[i+1,j+1] + pd*V[i+2,j+1]);
                    else
                    {
                        if(PutCall == "C")
                            V[i,j] = Math.Max(S[i,j] - K,Math.Exp(-r*dt)*(pu*V[i,j+1] + pm*V[i+1,j+1] + pd*V[i+2,j+1]));
                        else
                            V[i,j] = Math.Max(K - S[i,j],Math.Exp(-r*dt)*(pu*V[i,j+1] + pm*V[i+1,j+1] + pd*V[i+2,j+1]));
                    }
                }
            return V[0,0];
        }
    }
}

