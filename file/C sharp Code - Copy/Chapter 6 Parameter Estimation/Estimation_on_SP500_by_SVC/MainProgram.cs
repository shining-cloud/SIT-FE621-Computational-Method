using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Estimation_on_SP500_by_SVC
{
    class EstimationSVC
    {
        static void Main(string[] args)
        {
            // Classes
            HestonPrice HP = new HestonPrice();
            BlackScholesPrice BS = new BlackScholesPrice();
            Bisection BA = new Bisection();
            NelderMeadAlgo NM = new NelderMeadAlgo();
            ObjectiveFunction OF = new ObjectiveFunction();

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

            // Option settings
            OPSet opset;
            opset.S = 137.14;
            opset.r = 0.0010;
            opset.q = 0.0068;
            opset.trap = 1;

            // Read in SP500 implied volatilities
            int NT = 4;
            int NK = 7;
            double[,] MktIV = new Double[7,4] {
                {0.2780, 0.2638, 0.2532, 0.2518},{0.2477, 0.2402, 0.2364, 0.2369},
                {0.2186, 0.2158, 0.2203, 0.2239},{0.1878, 0.1930, 0.2047, 0.2098},
                {0.1572, 0.1712, 0.1894, 0.1970},{0.1334, 0.1517, 0.1748, 0.1849},
                {0.1323, 0.1373, 0.1618, 0.1736}};
            double[] K = new Double[7] { 120.0,125.0,130.0,135.0,140.0,145.0,150.0 };
            double[] T = new Double[4] {0.123287671232877, 0.268493150684932, 0.715068493150685, 0.953424657534247};

            // PutCall identifiers
            string[,] PutCall = new String[NK,NT];
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<=NT-1;t++)
                    PutCall[k,t] = "C";

            // Obtain the market prices
            double[,] MktPrice = new Double[NK,NT];
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<=NT-1;t++)
                    MktPrice[k,t] = BS.BlackScholes(opset.S,K[k],T[t],opset.r,opset.q,MktIV[k,t],PutCall[k,t]);

            // settings for the market data 
            MktData data = new MktData();
            data.MktIV   = MktIV;
            data.MktPrice  = MktPrice;
            data.K       = K;
            data.T       = T;
            data.PutCall = PutCall;

            // Bounds on the parameter estimates
            // kappa theta sigma v0 rho
            double e = 1e-5;
            double[] lb = new double[5] { e,e,e,e,-0.99 };
            double[] ub = new double[5] { 20.0,2.0,2.0,2.0,0.99 };

            // C.F. choice and Loss function choice
            string CF = "Heston";                   // Choice of c.f. "Heston" or "Attari"
            int LossFunction = 1;                   // Choice of loss function
                                                    // 1=MSE, 2=RMSE, 3=IVMSE, 4=Christoffersen et al.

            // Settings for the objective function
            OFSet ofset;
            ofset.opset = opset;
            ofset.data = data;
            ofset.X = X;
            ofset.W = W;
            ofset.LossFunction = LossFunction;
            ofset.lb = lb;
            ofset.ub = ub;
            ofset.CF = CF;

            // Settings for the Nelder Mead algorithm
            NMSet nmsettings;
            nmsettings.ofset = ofset;
            nmsettings.N = 5;						// Number of Heston parameters 
            nmsettings.MaxIters = 1000;		        // Maximum number of iterations
            nmsettings.Tolerance = 1e-4;		    // Tolerance on best and worst function values

            // Starting values (vertices) in vector form.  Add random increment about each starting value
            double kappaS = 9;
            double thetaS = 0.05;
            double sigmaS = 0.3;
            double v0S    = 0.05;
            double rhoS   = -0.8;
            int N = nmsettings.N;
            double[,] s = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                s[0,j] = kappaS + NM.RandomNum(-0.10,0.10);
                s[1,j] = thetaS + NM.RandomNum(-0.01,0.01);
                s[2,j] = sigmaS + NM.RandomNum(-0.05,0.05);
                s[3,j] = v0S    + NM.RandomNum(-0.01,0.01);
                s[4,j] = rhoS   + NM.RandomNum(-0.05,0.05);
            }

            // Find the Nelder-Mead parameter estimates
            double[] B = NM.NelderMead(OF.f,nmsettings,s);

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
            Console.WriteLine("Number of iterations required       {0}",B[6]);
            Console.WriteLine("  ");
            Console.WriteLine("----------------------------------------");

            // Obtain the model prices and model implied volatilities
            HParam paramNM = new HParam();
            paramNM.kappa = B[0];
            paramNM.theta = B[1];
            paramNM.sigma = B[2];
            paramNM.v0    = B[3];
            paramNM.rho   = B[4];
            double a = 0.01;
            double b = 3.0;
            double Tol = 1e-5;
            int MaxIter = 1000;

            double[,] ModelPrice = new double[NK,NT];
            double[,] ModelIV    = new double[NK,NT];
            double IVMSE = 0.0;

            double S = ofset.opset.S;
            double r = ofset.opset.r;
            double q = ofset.opset.q;
            int trap = ofset.opset.trap;

            for(int k=0;k<NK;k++)
                for(int t=0;t<NT;t++)
                {
                    ModelPrice[k,t] = HP.HestonPriceGaussLaguerre(paramNM,S,K[k],r,q,T[t],trap,PutCall[k,t],X,W);
                    ModelIV[k,t] = BA.BisecBSIV(PutCall[k,t],S,K[k],r,q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
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

    }
}

