using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Black_Scholes
{
    class MSBlackScholes
    {
        static void Main(string[] args)
        {
            // Reproduces Table 2 in Medvedev and Scaillet (2010) for American puts
            // under the Black Scholes model

            // Option settings
            double S = 40.0;
            double r = 0.0488;
            double q = 0.00;

            // Medvedev-Scaillet put option settings
            MSset mssettings = new MSset();
            mssettings.r = r;
            mssettings.q = q;

            // Trinomial tree settings
            int N = 500; 
            string PutCall  = "P";
            string EuroAmer = "A";

            // Table 2 settings
            double[] sigma = new double[9] { 0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4 };
            double[] T = new double[9] { 1.0/12.0,1.0/3.0,7.0/12.0,1.0/12.0,1.0/3.0,7.0/12.0,1.0/12.0,1.0/3.0,7.0/12.0 };
            double[,] BSPut = new double[3,9];     // Black-Scholes
            double[,] MSPut = new double[3,9];     // Medvedev-Scaillet
            double[,] BTPut = new double[3,9];     // Trinomial tree
            double SumError = 0.0;

            // Settings for the Bisection Method
            int MaxIter = 20000;
            double tol = 1e-20;
            double dt = 1e-10;
            double b =  3.0;

            // Find the tree and M-S prices
            BlackScholesPrice BS = new BlackScholesPrice();
            MSPrice MS = new MSPrice();
            TrinomialPrice TP = new TrinomialPrice();
            double[] K = new double[3] {35.0,40.0,45.0};
            for(int k=0;k<=2;k++)
            {
                for(int i=0;i<=8;i++)
                {
                    // Black Scholes European put
                    BSPut[k,i] = BS.BlackScholes(S,K[k],T[i],r,q,sigma[i],"P");
                    // Medvedev Scaillet American Put
                    MSPut[k,i] = MS.MSPriceBS(S,K[k],T[i],sigma[i],r,q,MaxIter,tol,b,dt);
                    // Trinomial Tree American put
                    BTPut[k,i] = TP.TrinomialTree(S,K[k],r,q,T[i],sigma[i],N,PutCall,EuroAmer);
                    SumError += Math.Abs(MSPut[k,i] - BTPut[k,i]);
                }
            }

            // Output the results for K = 35, 40, 45
            Console.WriteLine(" ");
            Console.WriteLine("                       Table 2 of Medvedev and Scaillet (2010)");
            Console.WriteLine("----------------------------------------------------------------------------");
            Console.WriteLine("                sigma = 0.2          sigma = 0.3            sigma = 0.4");
            Console.WriteLine("         ---------------------  ---------------------  ---------------------");
            Console.WriteLine("K = 35   T=1/12  T=1/3  T=7/12  T=1/12  T=1/3  T=7/12  T=1/12  T=1/3  T=7/12 ");
            Console.WriteLine("----------------------------------------------------------------------------");
            for(int k=0;k<=2;k++)
            {
                Console.WriteLine("EuroPut  {0,6:F3} {1,6:F3} {2,6:F3} {3,8:F3} {4,6:F3} {5,6:F3} {6,8:F3} {7,6:F3} {8,7:F3}",
                    BSPut[k,0],BSPut[k,1],BSPut[k,2],BSPut[k,3],BSPut[k,4],BSPut[k,5],BSPut[k,6],BSPut[k,7],BSPut[k,8]);
                Console.WriteLine("MSPut    {0,6:F3} {1,6:F3} {2,6:F3} {3,8:F3} {4,6:F3} {5,6:F3} {6,8:F3} {7,6:F3} {8,7:F3}",
                    MSPut[k,0],MSPut[k,1],MSPut[k,2],MSPut[k,3],MSPut[k,4],MSPut[k,5],MSPut[k,6],MSPut[k,7],MSPut[k,8]);
                Console.WriteLine("TreePut  {0,6:F3} {1,6:F3} {2,6:F3} {3,8:F3} {4,6:F3} {5,6:F3} {6,8:F3} {7,6:F3} {8,7:F3}",
                    BTPut[k,0],BTPut[k,1],BTPut[k,2],BTPut[k,3],BTPut[k,4],BTPut[k,5],BTPut[k,6],BTPut[k,7],BTPut[k,8]);
               Console.WriteLine("----------------------------------------------------------------------------");
            }
            Console.WriteLine("EuroPut = Black-Scholes closed form European put price");
            Console.WriteLine("MSPut   = Fifth-order MS (2001) American put expansion");
            Console.WriteLine("TreePut = Trinomial tree American put with {0} steps",N);
            Console.WriteLine("Sum of absolute errors {0,8:F5}",SumError);
            Console.WriteLine("----------------------------------------------------------------------------");
            Console.WriteLine(" ");
        }
    }
}
