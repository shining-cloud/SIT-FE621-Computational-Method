using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Differential_Evolution
{
    class DE
    {
        static void Main(string[] args)
        {
            // Classes
            MiscFunctions MF = new MiscFunctions();
            DiffEvoAlgo DE = new DiffEvoAlgo();
            NelderMeadAlgo NM = new NelderMeadAlgo();
            ObjectiveFunction OF = new ObjectiveFunction();

            // 32-point Gauss-Laguerre Abscissas and weights
            double[] X = new Double[32];
            double[] W = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    X[k] = double.Parse(bits[0]);
                    W[k] = double.Parse(bits[1]);
                }

            // Option settings
            OPSet opset;
            opset.S = 137.14;
            opset.r = 0.0010;
            opset.q = 0.0068;
            opset.trap = 1;

            // Bisection algorithm settings
            double a = 0.01;
            double b = 3.0;
            double Tol = 1e-5;
            int MaxIter = 1000;
            int ObjectiveFun = 1;

            // Read in SP500 implied volatilities
            int NT = 4;
            int NK = 7;
            double[,] MktIV = new Double[7,4] 
            {{0.2780, 0.2638, 0.2532 ,0.2518},
             {0.2477, 0.2402, 0.2364, 0.2369},
             {0.2186, 0.2158, 0.2203, 0.2239},
             {0.1878, 0.1930, 0.2047, 0.2098},
             {0.1572, 0.1712, 0.1894, 0.1970},
             {0.1334, 0.1517, 0.1748, 0.1849},
             {0.1323, 0.1373, 0.1618, 0.1736}};
            double[] K = new Double[7] { 120.0,125.0,130.0,135.0,140.0,145.0,150.0 };
            double[] T = new Double[4] { 0.123287671232877,0.268493150684932,0.715068493150685,0.953424657534247 };

            // PutCall identifiers
            string[,] PutCall = new String[NK,NT];
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<=NT-1;t++)
                    PutCall[k,t] = "C";

            // Obtain the market prices
            BlackScholesPrice BS = new BlackScholesPrice();
            double[,] MktPrice = new Double[NK,NT];
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<=NT-1;t++)
                    MktPrice[k,t] = BS.BlackScholes(opset.S,K[k],T[t],opset.r,opset.q,MktIV[k,t],PutCall[k,t]);

            // Place the market data in the structure
            MktData data = new MktData();
            data.MktIV   = MktIV;
            data.MktPrice  = MktPrice;
            data.K       = K;
            data.T       = T;
            data.PutCall = PutCall;

            // Estimation bounds
            double e = 1e-5;
            double kappaL = e; double kappaU = 10;
            double thetaL = e; double thetaU = 5;
            double sigmaL = e; double sigmaU = 5;
            double v0L    = e; double v0U    = 1;
            double rhoL   = -.9; double rhoU   = 0;
            double[] ub = new double[5] { kappaU,thetaU,sigmaU,v0U,rhoU };
            double[] lb = new double[5] { kappaL,thetaL,sigmaL,v0L,rhoL };

            // Objective function settings;
            OFSet ofset;
            ofset.opset = opset;
            ofset.data = data;
            ofset.X = X;
            ofset.W = W;
            ofset.LossFunction = 4;           // Choice of loss function 1=MSE, 2=RMSE, 3=IVMSE, 4=Christoffersen et al.
            ofset.lb = lb;
            ofset.ub = ub;
            ofset.CF = "Heston";           // Choice of c.f. "Heston" or "Attari"

            // Settings for the Differential Evolution algorithm
            DEParam ParamLim = new DEParam();
            ParamLim.ub = ub;
            ParamLim.lb = lb;
            ParamLim.NG = 500;
            ParamLim.NP = 75;
            ParamLim.F  = 0.5;
            ParamLim.CR = 0.8;

            // Run the differential evolution algorithm
            HParam DEparam = DE.HestonDE(ParamLim,opset,data,ObjectiveFun,a,b,Tol,MaxIter,X,W,ofset.CF);

            // Starting values (vertices) in vector form.  Add random increment about each starting value
            double kappaS = 9;
            double thetaS = 0.05;
            double sigmaS = 0.3;
            double v0S    = 0.05;
            double rhoS   = -0.8;
            int N = 5;
            double[,] s = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                s[0,j] = kappaS + MF.RandomNum(-0.10,0.10);
                s[1,j] = thetaS + MF.RandomNum(-0.01,0.01);
                s[2,j] = sigmaS + MF.RandomNum(-0.05,0.05);
                s[3,j] = v0S    + MF.RandomNum(-0.01,0.01);
                s[4,j] = rhoS   + MF.RandomNum(-0.05,0.05);
            }

            // Nelder Mead settings
            NMSet nmset;
            nmset.ofset = ofset;
            nmset.N = 5;						// Number of Heston parameters 
            nmset.MaxIters = 1000;		    // Maximum number of iterations
            nmset.Tolerance = 1e-4;		// Tolerance on best and worst function values

            // Run the Nelder Mead algorithm
            double[] B = NM.NelderMead(OF.f,nmset,s);
            HParam NMparam = new HParam();
            NMparam.kappa = B[0];
            NMparam.theta = B[1];
            NMparam.sigma = B[2];
            NMparam.v0    = B[3];
            NMparam.rho   = B[4];

            // Calculate IVMSE under both parameter estimates
            double[,] NMPrice = new double[NK,NT];
            double[,] DEPrice = new double[NK,NT];
            double[,] NMIV = new double[NK,NT];
            double[,] DEIV = new double[NK,NT];
            double NMIVMSE = 0.0;
            double DEIVMSE = 0.0;
            double S = opset.S;
            double r = opset.r;
            double q = opset.q;
            int trap = opset.trap;

            HestonPrice HP = new HestonPrice();
            Bisection BA = new Bisection();
            for(int t=0;t<NT;t++)
            {
                for(int k=0;k<NK;k++)
                {
                    NMPrice[k,t] = HP.HestonPriceGaussLaguerre(NMparam,S,K[k],r,q,T[t],trap,PutCall[k,t],X,W);
                    NMIV[k,t] = BA.BisecBSIV(PutCall[k,t],S,K[k],r,q,T[t],a,b,NMPrice[k,t],Tol,MaxIter);
                    NMIVMSE += Math.Pow(MktIV[k,t] - NMIV[k,t],2) / Convert.ToDouble(NT*NK);
                    DEPrice[k,t] = HP.HestonPriceGaussLaguerre(DEparam,S,K[k],r,q,T[t],trap,PutCall[k,t],X,W);
                    DEIV[k,t] = BA.BisecBSIV(PutCall[k,t],S,K[k],r,q,T[t],a,b,DEPrice[k,t],Tol,MaxIter);
                    DEIVMSE += Math.Pow(MktIV[k,t] - DEIV[k,t],2) / Convert.ToDouble(NT*NK);
                }
            }

            // Output the result
            Console.WriteLine("  ");
            Console.WriteLine("Parameter Estimates --------------------");
            Console.WriteLine("                 Differential Evolution     Nelder Mead");
            Console.WriteLine("kappa            {0,10:F4} {1,25:F4}",DEparam.kappa,NMparam.kappa);
            Console.WriteLine("theta            {0,10:F4} {1,25:F4}",DEparam.theta,NMparam.theta);
            Console.WriteLine("sigma            {0,10:F4} {1,25:F4}",DEparam.sigma,NMparam.sigma);
            Console.WriteLine("v0               {0,10:F4} {1,25:F4}",DEparam.v0,NMparam.v0);
            Console.WriteLine("rho              {0,10:F4} {1,25:F4}",DEparam.rho,NMparam.rho);
            Console.WriteLine("  ");
            Console.WriteLine("IV MSE ---------------------------------");
            Console.WriteLine("Differential Evolution IVMSE {0:E8}",DEIVMSE);
            Console.WriteLine("Nelder Mead IVMSE            {0:E8}",NMIVMSE);
            Console.WriteLine("  ");
        }
    }
}
