using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_LSM
{
    class LSMonteCarlo
    {
        // Longstaff-Schwartz for American puts and calls
        // Also returns Euro price
        public double[] HestonLSM(double[,] S,double K,double r,double q,double T,int NT,int NS,string PutCall)
        {
            Regression RE = new Regression();

            // Time increment
            double dt = T/Convert.ToDouble(NT);

            // Initialize the Cash Flows.
            double[,] CF = new double[NS,NT];

            // Set the last cash flows to the intrinsic value.
            for(int s=0;s<=NS-1;s++)
                if(PutCall == "P")
                    CF[s,NT-1] = Math.Max(K - S[s,NT-1],0.0);
                else if(PutCall == "C")
                    CF[s,NT-1] = Math.Max(S[s,NT-1] - K,0.0);

            // European price
            double EuroPrice = 0.0;
            for(int s=0;s<=NS-1;s++)
                EuroPrice += Math.Exp(-r*T)*CF[s,NT-1]/Convert.ToDouble(NS);

            // Work backwards through the stock prices until time t=2.
            // We could work through to time t=1 but the regression will not be
            // of full rank at time 1, so this is cleaner.
            for(int t=NT-2;t>=1;t--)
            {
                // Indices for stock paths in-the-money at time t
                int[] I = new int[NS];
                for(int s=0;s<=NS-1;s++)
                {
                    I[s] = 0;
                    if(((PutCall == "P") & (S[s,t] < K)) | ((PutCall == "C") & (S[s,t] > K)))
                        I[s] = 1;
                }

                // Stock paths in-the-money at time t
                int NI = 0;
                List<double> X  = new List<double>();
                List<int> Xi = new List<int>();
                for(int s=0;s<=NS-1;s++)
                    if(I[s] == 1)
                    {
                        X.Add(S[s,t]);
                        Xi.Add(s);
                        NI += 1;
                    }

                // Cash flows at time t+1, discounted one period
                double[] Y = new double[NI];
                for(int s=0;s<=NI-1;s++)
                    Y[s] = CF[Xi[s],t+1]*Math.Exp(-r*dt);

                // Design matrix for regression to predict cash flows
                double[,] Z = new double[NI,3];

                for(int s=0;s<=NI-1;s++)
                {
                    Z[s,0] = 1.0;
                    Z[s,1] = (1.0 - X[s]);
                    Z[s,2] = (2.0 - 4.0*X[s] - X[s]*X[s])/2.0;
                }

                // Regression parameters and predicted cash flows
                double[] beta = RE.Beta(Z,Y);
                double[] PredCF = RE.MVMult(Z,beta);

                // Indices for stock paths where immediate exercise is optimal
                // J[s] contains the path number
                // E[s] contains the 
                List<int> E = new List<int>();
                List<int> J = new List<int>();
                int NE = 0;
                for(int s=0;s<=NI-1;s++)
                    if((PutCall == "P" & (K - X[s]>0  &  K - X[s]> PredCF[s])) |
                       (PutCall == "C" & (X[s] - K>0  &  X[s] - K> PredCF[s])))
                    {
                        J.Add(s);
                        E.Add(Xi[s]);
                        NE += 1;
                    }

                // All other paths --> Continuation is optimal
                int[] All  = new int[NS];
                for(int k=0;k<=NS-1;k++) All[k] = k;

                // C contains indices for continuation paths
                IEnumerable<int> C = All.Except(E);
                int[] Co = C.ToArray();
                int NC = Co.Length;

                // Replace cash flows with exercise value where exercise is optimal
                for(int s=0;s<=NE-1;s++)
                    if(PutCall == "P")
                        CF[E[s],t] = Math.Max(K - X[J[s]],0.0);
                    else if(PutCall == "C")
                        CF[E[s],t] = Math.Max(X[J[s]] - K,0.0);
                for(int s=0;s<=NC-1;s++)
                    CF[Co[s],t] = Math.Exp(-r*dt)*CF[Co[s],t+1];
            }
            // Calculate the cash flow at time 2
            double[] Op = new double[NS];
            for(int s=0;s<=NS-1;s++)
                Op[s] = CF[s,1];

            // Calculate the American price
            double AmerPrice = Math.Exp(-r*dt)*RE.VMean(Op);
            
            // Return the European and American prices
            double[] output = new double[2] { EuroPrice,AmerPrice };
            return output;
        }
    }
}

