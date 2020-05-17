using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Beliaeva_Nawalkha_Tree
{
    class VolStruct:IDisposable
    {
        public float[,] X;
        public float[,] V;
        public int RBound;
        public int M;
        bool disposed;

        protected virtual void Dispose(bool disposing)
        {
            if(!disposed)
            {
                if(disposing)
                {
                    //dispose managed ressources
                }
            }
            //dispose unmanaged ressources
            disposed = true;
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
    }

    class VolTree
    {
        public VolStruct BuildVolTree(float kappa,float theta,float sigma,float V0,float dt,int NT,float threshold)
        {
            // Creates the Beliaeva-Nawalkha tree for the variance process
            // INPUTS
            //   kappa = Heston parameter for mean reversion speed
            //   theta = Heston parameter for mean reversion level
            //   sigma = Heston parameter for vol of vol
            //   V0     = Heston parameter for initial variance
            //   rf = Risk free rate
            //   dt = time increment
            //   NT = Number of time steps
            //   threshold = threshold for which V can be zero
            // OUTPUTS
            //   X  = the tree for transformed variance
            //   V  = the tree for the original variance
            //   RBound = the row where the V(n,t) variances are zero
            //   M = row where volatility starts at time 1

            // Initial quantities
            float X0 = Convert.ToSingle(2.0*Math.Sqrt(V0)/sigma);

            float be = Convert.ToSingle(X0/Math.Sqrt(dt)/Math.Floor(X0/Math.Sqrt(1.5*dt)));
            float bc = Convert.ToSingle(X0/Math.Sqrt(dt)/Math.Floor(X0/Math.Sqrt(1.5*dt)+1.0));
            float b;
            if(Math.Abs(bc-Math.Sqrt(1.5)) < Math.Abs(be-Math.Sqrt(1.5)))
                b = bc;
            else
                b = be;
            if((b<1.0) || (b>Math.Sqrt(2.0)))
                Console.WriteLine("Warning b = {0:F5} but should be 1 < b < 1.4142\n",b);

            // Initialize the X-matrix and V-matrix
            int NR = 2*NT - 1;
            float[,] X = new float[NR,NT];
            float[,] V = new float[NR,NT];
            int M = (NR+1)/2 - 1;
            X[M,0] = X0;

            // Time 1 node for X
            float muX = Convert.ToSingle(1.0/X0*(0.5*kappa*(4.0*theta/sigma/sigma-X0*X0)-0.5));
            float J = Convert.ToSingle(Math.Floor(muX*Math.Sqrt(dt)/b + 1.0/b/b));
            X[M-1,1] = Convert.ToSingle(X[M,0] + b*(J+1.0)*Math.Sqrt(dt));
            X[M+0,1] = Convert.ToSingle(X[M,0] + b*(J+0.0)*Math.Sqrt(dt));
            X[M+1,1] = Convert.ToSingle(X[M,0] + b*(J-1.0)*Math.Sqrt(dt));

            // Remaining nodes for X
            for(int t=1;t<=NT-2;t++)
            {
                for(int n=M-t;n<=M+t;n++)
                {
                    muX = Convert.ToSingle(1.0/X[n,t]*(0.5*kappa*(4.0*theta/sigma/sigma-X[n,t]*X[n,t])-0.5));
                    J = Convert.ToSingle(Math.Floor(muX*Math.Sqrt(dt)/b + 1/b/b));
                    // Nodes where X > 0
                    if(X[n,t]>threshold & X[n,t]*X[n,t]*sigma*sigma/4.0>threshold)
                    {
                        X[n-1,t+1] = Convert.ToSingle(X[n,t] + b*(J+1.0)*Math.Sqrt(dt));
                        X[n+0,t+1] = Convert.ToSingle(X[n,t] + b*(J+0.0)*Math.Sqrt(dt));
                        X[n+1,t+1] = Convert.ToSingle(X[n,t] + b*(J-1.0)*Math.Sqrt(dt));
                    }
                    // Nodes where X = 0
                    else
                    {
                        X[n-1,t+1] = X[n-1,t];
                        X[n+0,t+1] = X[n+0,t];
                        X[n+1,t+1] = X[n+1,t];
                    }
                }
            }

            // Identify row of the tree where X = 0
            int RBound = 0;
            while(X[RBound,M]>=threshold && RBound < NR-1)
                RBound += 1;

            //// Build the volatility tree V(n,t)
            for(int t=0;t<=NT-1;t++)
                for(int n=0;n<=NR-1;n++)
                {
                    V[n,t] = Convert.ToSingle(X[n,t]*X[n,t]*sigma*sigma/4.0);
                    if(V[n,t] < threshold)
                        V[n,t] = 0.0f;
                    if(X[n,t] < threshold)
                        X[n,t] = 0.0f;
                }

            // Output the results
            VolStruct output = new VolStruct();
            output.X = X;
            output.V = V;
            output.M = M;
            output.RBound = RBound;
            return output;
        }
    }
}


