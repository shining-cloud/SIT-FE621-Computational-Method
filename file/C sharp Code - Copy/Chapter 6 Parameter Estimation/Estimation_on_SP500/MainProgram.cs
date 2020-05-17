using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Estimation_on_SP500
{
    class EstimationSP500
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

            // Bounds on the parameter estimates
            // kappa theta sigma v0 rho
            double e = 1e-5;
            double[] lb = new double[5] { e,e,e,e,-0.99 };
            double[] ub = new double[5] { 20.0,2.0,2.0,2.0,0.99 };

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
            double[] T = new Double[4] {0.123287671232877, 0.268493150684932, 0.715068493150685, 0.953424657534247};

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

            // Settings for the objective function
            OFSet ofsettings;
            ofsettings.opsettings = opsettings;
            ofsettings.data = data;
            ofsettings.X = X;
            ofsettings.W = W;
            ofsettings.LossFunction = 1;
            ofsettings.lb = lb;
            ofsettings.ub = ub;

            // Settings for the Nelder Mead algorithm
            NMSet nmsettings;
            nmsettings.N = 5;						// Number of Heston parameters 
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
            double[,] x = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                x[0,j] = kappaS + RandomNum(-0.01,0.01);
                x[1,j] = thetaS + RandomNum(-0.01,0.01);
                x[2,j] = sigmaS + RandomNum(-0.01,0.01);
                x[3,j] = v0S    + RandomNum(-0.01,0.01);
                x[4,j] = rhoS   + RandomNum(-0.01,0.01);
            }
            
            // Obtain the parameter estimates
            NelderMeadAlgo NM = new NelderMeadAlgo();
            ObjectiveFunction OF = new ObjectiveFunction();
            double[] B = NM.NelderMead(OF.f,nmsettings,x);

            // Output the estimation result
            Console.WriteLine("  ");
            Console.WriteLine("Parameter Estimates --------------------");
            Console.WriteLine("  ");
            Console.WriteLine("kappa   =  {0:F5}",B[0]);
            Console.WriteLine("theta   =  {0:F5}",B[1]);
            Console.WriteLine("sigma   =  {0:F5}",B[2]);
            Console.WriteLine("v0      =  {0:F5}",B[3]);
            Console.WriteLine("rho     =  {0:F5}",B[4]);
            Console.WriteLine("  ");
            Console.WriteLine("Value of the objective function is  {0:F5}",B[5]);
            Console.WriteLine("  ");
            Console.WriteLine("Number of iterations required       {0:0}",B[6]);
            Console.WriteLine("  ");
            Console.WriteLine("----------------------------------------");

            // Obtain the model prices and model implied volatilities
            HParam paramNM = new HParam();
            paramNM.kappa = B[0];
            paramNM.theta = B[1];
            paramNM.sigma = B[2];
            paramNM.v0    = B[3];
            paramNM.rho   = B[4];

            // Settings for bisection algorithm
            double a = 0.01;
            double b = 3.00;
            double Tol = 1e-4;
            int MaxIter = 1000;

            double[,] ModelPrice = new double[NK,NT];
            double[,] ModelIV    = new double[NK,NT];
            double IVMSE = 0.0;

            HestonPrice HP = new HestonPrice();
            Bisection BA = new Bisection();
            for(int k=0;k<NK;k++)
                for(int t=0;t<NT;t++)
                {
                    ModelPrice[k,t] = HP.HestonPriceGaussLaguerre(paramNM,opsettings.S,K[k],opsettings.r,opsettings.q,T[t],opsettings.trap,PutCall[k,t],X,W);
                    ModelIV[k,t] = BA.BisecBSIV(PutCall[k,t],opsettings.S,K[k],opsettings.r,opsettings.q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                    IVMSE += Math.Pow(MktIV[k,t] - ModelIV[k,t],2) / Convert.ToDouble(NT*NK);
                }

            // Output the results
            Console.Write("MSE between model and market implied vols  = {0}",IVMSE);
            Console.WriteLine("  ");
            Console.WriteLine("----------------------------------------");
            Console.WriteLine("Market implied volatilities");
            Console.WriteLine("  ");
            for(int k=0;k<NK;k++)
                for(int t=0;t<NT;t++)
                    if(t<NT-1)
                        Console.Write("{0:F4}   ",MktIV[k,t]);
                    else
                        Console.WriteLine("{0:F4}   ",MktIV[k,t]);

            Console.WriteLine("  ");
            Console.WriteLine("----------------------------------------");
            Console.WriteLine("Model implied volatilities");
            Console.WriteLine("  ");
            for(int k=0;k<NK;k++)
                for(int t=0;t<NT;t++)
                    if(t<NT-1)
                        Console.Write("{0:F4}   ",ModelIV[k,t]);
                    else
                        Console.WriteLine("{0:F4}   ",ModelIV[k,t]);

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

