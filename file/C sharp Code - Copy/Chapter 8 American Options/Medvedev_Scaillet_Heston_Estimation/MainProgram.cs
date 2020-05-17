using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Medvedev_Scaillet_Heston
{
    class MedvedevScailletEstimation
    {
        static void Main(string[] args)
        {
            // 32 point Gauss Laguerre
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
            double[] ub = new double[5] { 20.0,2.0,5.0,2.0,0.99 };

            //// IBM put prices May 7, 2010
            double S = 122.10;
            //double[,] MktPrice = new double[8,5]{
            //    { 0.130,  0.620,  1.275,  2.950,  4.750},
            //    { 0.260,  0.955,  1.830,  3.925,  6.075},
            //    { 0.485,  1.500,  2.610,  5.200,  7.625},
            //    { 0.995,  2.445,  3.775,  6.800,  9.375},
            //    { 2.155,  3.900,  5.475,  8.800, 11.550},
            //    { 4.525,  6.225,  7.775, 11.225, 13.975},
            //    { 8.375,  9.525, 10.850, 14.125, 16.775},
            //    {13.075, 13.600, 14.575, 17.425, 19.900}};

            // Use the first set of prices, for the first maturity
            double[] MktPrice = new double[8] {0.130,0.260,0.485,0.995,2.155,4.525,8.375,13.075};
            double[] K = new double[8] {100.0,105.0,110.0,115.0,120.0,125.0,130.0,135.0};
            double r = 0.015;
            double q = 0.01;
            //double[] T = new double[5] {0.0384,0.1151, 0.1918, 0.4411, 0.7096};
            double T = 0.0384;
            int NK = 8;

            // Bisection algorithm
            double a = 0.01;
            double b = 2.00;
            double Tol = 1.0e-5;
            int MaxIter = 1000;
            double B = 2.0;
            double dt = 1.0e-10;

            // Implied volatilities
            BisectionImpliedVol BA = new BisectionImpliedVol();
            double[] MktIV = new double[NK];
            for(int k=0;k<=NK-1;k++)
            {
                MktIV[k] = BA.BisectionMSIV(S,K[k],r,q,T,a,b,MktPrice[k],Tol,MaxIter,B,dt);
            }

            MktData MktData = new MktData();
            MktData.MktIV = MktIV;
            MktData.K     = K;
            MktData.PutCall = "P";

            OpSet opsettings = new OpSet();
            opsettings.S = S;
            opsettings.T = T;
            opsettings.r = r;
            opsettings.q = q;
            opsettings.PutCall = "P";
            opsettings.trap = 1;

            MSSet mssettings = new MSSet();
            mssettings.method = 3;
            mssettings.A = 0.0001;
            mssettings.B = 100.0;
            mssettings.N = 3000;
            mssettings.dt  = 1.0e-10;
            mssettings.tol = 1.0e-5;
            mssettings.MaxIter = 1000;
            mssettings.NumTerms = 3;
            mssettings.yinf = 1.0e4;
            mssettings.a = -10.0;
            mssettings.b = +10.0;

            OFSet ofsettings = new OFSet();
            ofsettings.opsettings = opsettings;
            ofsettings.mssettings = mssettings;
            ofsettings.data       = MktData;
            ofsettings.X = X;
            ofsettings.W = W;
            ofsettings.lb = lb;
            ofsettings.ub = ub;

            // Settings for the Nelder Mead algorithm
            NMSet nmsettings;
            nmsettings.N = 5;					// Number of Heston parameters 
            nmsettings.MaxIters = 20;		     // Maximum number of iterations
            nmsettings.Tolerance = 1e-6;		// Tolerance on best and worst function values
            nmsettings.ofsettings = ofsettings;

            // Starting values (vertices) in vector form.  Add random increment about each starting value
            NelderMeadAlgo NM = new NelderMeadAlgo();
            double kappaS =  20.00;
            double thetaS =  0.036;
            double sigmaS =  3.9;
            double v0S    =  0.15;
            double rhoS   = -0.46;
            int N = nmsettings.N;
            double[,] x = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                x[0,j] = kappaS + NM.RandomNum(-0.01,0.01)*kappaS;
                x[1,j] = thetaS + NM.RandomNum(-0.01,0.01)*thetaS;
                x[2,j] = sigmaS + NM.RandomNum(-0.01,0.01)*sigmaS;
                x[3,j] = v0S    + NM.RandomNum(-0.01,0.01)*v0S;
                x[4,j] = rhoS   + NM.RandomNum(-0.01,0.01)*rhoS;
            }

            // Obtain the parameter estimates
            ObjectiveFunction OF = new ObjectiveFunction();
            double[] BB = NM.NelderMead(OF.f,nmsettings,x);

            HParam paramEst = new HParam();
            paramEst.kappa = BB[0];
            paramEst.theta = BB[1];
            paramEst.sigma = BB[2];
            paramEst.v0    = BB[3];
            paramEst.rho   = BB[4];
            paramEst.lambda = 0.0;

            // Obtain the fitted implied vols
            MSPrices MS = new MSPrices();
            double[] ModelIV    = new double[NK];
            double[] ModelPrice = new double[NK];
            double[] output     = new double[6];
            double IVMSE = 0.0;
            for(int k=0;k<=NK-1;k++)
            {
                opsettings.K = K[k];
                output = MS.MSPriceHeston(paramEst,opsettings,mssettings);
                ModelPrice[k] = output[2];
                ModelIV[k]    = BA.BisectionMSIV(S,K[k],r,q,T,a,b,ModelPrice[k],Tol,MaxIter,B,dt);
                IVMSE += Math.Pow(ModelIV[k] - MktIV[k],2.0);
            }

            // Output the estimation result
            Console.WriteLine("  ");
            Console.WriteLine("Parameter Estimates --------------------");
            Console.WriteLine("  ");
            Console.WriteLine("kappa   =  {0:F5}",BB[0]);
            Console.WriteLine("theta   =  {0:F5}",BB[1]);
            Console.WriteLine("sigma   =  {0:F5}",BB[2]);
            Console.WriteLine("v0      =  {0:F5}",BB[3]);
            Console.WriteLine("rho     =  {0:F5}",BB[4]);
            Console.WriteLine("  ");
            Console.WriteLine("Value of the objective function is  {0:E5}",BB[5]);
            Console.WriteLine("Number of iterations required       {0:0}",BB[6]);
            Console.WriteLine("  ");
            Console.WriteLine("----------------------------------------");
            Console.WriteLine(" ");
            Console.WriteLine("Strike    MktIV   ModelIV");
            Console.WriteLine("-------------------------");
            for(int k=0;k<=NK-1;k++)
                Console.WriteLine("{0,3:F0} {1,12:F5} {2,12:F5}",K[k],MktIV[k],ModelIV[k]);
        }
    }
}

