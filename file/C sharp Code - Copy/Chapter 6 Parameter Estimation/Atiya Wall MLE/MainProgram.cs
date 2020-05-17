using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;


namespace Atiya_Wall_MLE
{
    class AWMLE
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] X = new Double[32];
            double[] W = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    X[k] = double.Parse(bits[0]);
                    W[k] = double.Parse(bits[1]);
                }
            }
            // SP500 Historical prices.  Oldest prices are first
            double[] S = new Double[120];
            double[] x = new double[120];
            using(TextReader reader = File.OpenText("../../SPY Historical Prices.txt"))
            {
                for(int k=0;k<=119;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    S[k] = double.Parse(bits[0]);
                    x[k] = Math.Log(S[k]);
                }
            }

            // Bounds on the parameter estimates
            // kappa theta sigma v0 rho
            double e = 1e-2;
            double[] lb = new double[5] { e,e,e,e,-0.999 };
            double[] ub = new double[5] { 30.0,3.0,3.0,2.0,0.999 };

            // Option settings
            OPSet opsettings;
            opsettings.S = 137.14;
            opsettings.r = 0.0010;
            opsettings.q = 0.0068;
            opsettings.trap = 1;

            // Read in SP500 implied volatilities
            int NT = 4;
            int NK = 7;
            double[,] MktIV = new Double[7,4] {
                {0.2780, 0.2638, 0.2532, 0.2518}, {0.2477, 0.2402, 0.2364, 0.2369},
                {0.2186, 0.2158, 0.2203, 0.2239}, {0.1878, 0.1930, 0.2047, 0.2098},
                {0.1572, 0.1712, 0.1894, 0.1970}, {0.1334, 0.1517, 0.1748, 0.1849},
                {0.1323, 0.1373, 0.1618, 0.1736}};
            double[] K = new Double[7] { 120.0,125.0,130.0,135.0,140.0,145.0,150.0 };
            double[] T = new Double[4] { 45.0/365.0,98.0/365.0,261.0/365.0,348.0/365.0 };

            // PutCall identifiers
            string[,] PutCall = new String[NK,NT];
            for(int k=0;k<=NK-1;k++)
            {
                for(int t=0;t<=NT-1;t++)
                    PutCall[k,t] = "C";
            }

            // Obtain the market prices
            BlackScholesPrice BS = new BlackScholesPrice();
            double[,] MktPrice = new Double[NK,NT];
            for(int k=0;k<=NK-1;k++)
            {
                for(int t=0;t<=NT-1;t++)
                    MktPrice[k,t] = BS.BlackScholes(opsettings.S,K[k],T[t],opsettings.r,opsettings.q,MktIV[k,t],PutCall[k,t]);
            }

            // Place the market data in the structure
            MktData data = new MktData();
            data.MktIV   = MktIV;
            data.MktPrice  = MktPrice;
            data.K       = K;
            data.T       = T;
            data.PutCall = PutCall;

            // True parameter values
            int Lmethod = 2;
            double[] True = new double[5] { 8.9832,0.0524,1.0982,0.0325,-0.9921 };
            double dt = 1.0/252.0;

            // Settings for the objective function
            OFSet ofsettings;
            ofsettings.x = x;
            ofsettings.r = opsettings.r;
            ofsettings.q = opsettings.q;
            ofsettings.dt = dt;
            ofsettings.method = Lmethod;

            // Settings for the Nelder Mead algorithm
            NMSet nmsettings;
            nmsettings.N = 5;				    // Number of Heston parameters 
            nmsettings.MaxIters = 1000;		    // Maximum number of iterations
            nmsettings.Tolerance = 1e-3;		// Tolerance on best and worst function values
            nmsettings.ofsettings = ofsettings;

            // Starting values (vertices) in vector form.  Add random increment about each starting value
            double kappaS =  9.00;
            double thetaS =  0.05;
            double sigmaS =  0.30;
            double v0S    =  0.05;
            double rhoS   = -0.80;
            int N = nmsettings.N;
            double[,] xs = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                xs[0,j] = kappaS + RandomNum(-0.01,0.01);
                xs[1,j] = thetaS + RandomNum(-0.01,0.01);
                xs[2,j] = sigmaS + RandomNum(-0.01,0.01);
                xs[3,j] = v0S    + RandomNum(-0.01,0.01);
                xs[4,j] = rhoS   + RandomNum(-0.01,0.01);
            }

            // Obtain the parameter estimates
            NelderMeadAlgo NM = new NelderMeadAlgo();
            Likelihood LL = new Likelihood();
            double[] B = NM.NelderMead(LL.f,nmsettings,xs);

            // Output the estimation result
            Console.WriteLine("  ");
            Console.WriteLine("Atiya-Wall (2009) MLE parameters --------------------");
            Console.WriteLine("  ");
            Console.WriteLine("Parameter    MLE     True Value  ");
            Console.WriteLine("----------------------------------------");
            Console.WriteLine("kappa    {0,10:F5} {1,10:F5}",B[0],True[0]);
            Console.WriteLine("theta    {0,10:F5} {1,10:F5}",B[1],True[1]);
            Console.WriteLine("sigma    {0,10:F5} {1,10:F5}",B[2],True[2]);
            Console.WriteLine("v0       {0,10:F5} {1,10:F5}",B[3],True[3]);
            Console.WriteLine("rho      {0,10:F5} {1,10:F5}",B[4],True[4]);
            Console.WriteLine("----------------------------------------");
            Console.WriteLine("  ");
            Console.WriteLine("Value of the objective function is  {0:F5}",B[5]);
            Console.WriteLine("  ");
            Console.WriteLine("Number of iterations required       {0:0}",B[6]);
            Console.WriteLine("  ");
            Console.WriteLine("----------------------------------------");
        }
        // Random number in (a,b) ==========================================================
        private static readonly Random U = new Random();
        private static readonly object sync = new object();
        public static double RandomNum(double a,double b)
        {
            int divisor = 1000000000;
            lock(sync) { return a + (b-a)*U.Next(0,divisor)/divisor; }
        }
    }
}



