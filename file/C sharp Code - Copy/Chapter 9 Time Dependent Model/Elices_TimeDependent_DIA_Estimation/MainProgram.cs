using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Elices_TimeDependent_DIA_Estimation
{
    class ElicesEstimation
    {
        static void Main(string[] args)
        {
            // Classes
            BisectionAlgo BA = new BisectionAlgo();
            BlackScholesPrice BS = new BlackScholesPrice();
            HestonPrice HP = new HestonPrice();
            ElicesAlgo EA = new ElicesAlgo();
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

            // Bounds on the parameter estimates
            // kappa theta sigma v0 rho
            double e = 1.0e-2;
            double[] lb = new double[5] { e,e,e,-0.99,e };
            double[] ub = new double[5] { 10.0,3.0,10.0,0.99,3.0 };

            // Read in SP500 implied volatilities
            const int NT = 4;
            const int NK = 13;
            double[][] MktIV = new double[NT][];
            MktIV[0] = new double[NK] { 0.1962,0.1910,0.1860,0.1810,0.1761,0.1718,0.1671,0.1644,0.1645,0.1661,0.1701,0.1755,0.1796 };
            MktIV[1] = new double[NK] { 0.1947,0.1905,0.1861,0.1812,0.1774,0.1743,0.1706,0.1671,0.1641,0.1625,0.1602,0.1610,0.1657 };
            MktIV[2] = new double[NK] { 0.2019,0.1980,0.1943,0.1907,0.1871,0.1842,0.1813,0.1783,0.1760,0.1743,0.1726,0.1716,0.1724 };
            MktIV[3] = new double[NK] { 0.2115,0.2082,0.2057,0.2021,0.2000,0.1974,0.1950,0.1927,0.1899,0.1884,0.1862,0.1846,0.1842 };
            double[] K = new double[NK] { 124.0,125.0,126.0,127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0,136.0, };
            double[] T = new double[NT] { 37.0/365.0,72.0/365.0,135.0/365.0,226.0/365.0 };

            // PutCall identifier
            string PutCall = "P";

            // Maturity increments
            double[] tau = new double[NT];
            tau[0] = T[0];
            for(int t=1;t<=NT-1;t++)
                tau[t] = T[t] - T[t-1];

            // Option settings
            OPSet opsettings;
            opsettings.S = 129.14;
            opsettings.r = 0.0010;
            opsettings.q = 0.0068;
            opsettings.trap = 1;

            // Obtain the market prices for Elices, jagged array, and Heston (two dimensional array)
            double[][] MktPrice = new double[4][];
            double[,] MktPriceH = new double[NT,NK];
            double[,] MktIVH    = new double[NT,NK];
            for(int t=0;t<=NT-1;t++)
            {
                MktPrice[t] = new double[13];
                for(int k=0;k<=NK-1;k++)
                {
                    MktPrice[t][k] = BS.BlackScholes(opsettings.S,K[k],T[t],opsettings.r,opsettings.q,MktIV[t][k],PutCall);
                    MktPriceH[t,k] =  MktPrice[t][k];
                    MktIVH[t,k] = MktIV[t][k];
                }
            }

            // Settings for the Elices objective function
            OFSetE ofsetE = new OFSetE();
            ofsetE.opsettings = opsettings;
            ofsetE.X = X;
            ofsetE.W = W;
            ofsetE.LossFunction = 1;
            ofsetE.lb = lb;
            ofsetE.ub = ub;

            // Elices Settings for the Market data
            MktData MktData = new MktData();
            MktData.K = K;
            MktData.PutCall = PutCall;
            ofsetE.data = MktData;

            // Settings for the Heston objective function
            OFSetH ofsetH = new OFSetH();
            ofsetH.opsettings = opsettings;
            ofsetH.X = X;
            ofsetH.W = W;
            ofsetH.LossFunction = ofsetE.LossFunction;
            ofsetH.lb = ofsetE.lb;
            ofsetH.ub = ofsetE.ub;
            ofsetH.MktPrice = MktPriceH;
            ofsetH.MktIV    = MktIVH;
            ofsetH.K = K;
            ofsetH.T = T;
            ofsetH.PutCall = PutCall;

            // Parameter initial values
            double kappaS =  2.00;
            double thetaS =  0.10;
            double sigmaS =  1.20;
            double v0S    =  0.03;
            double rhoS   = -0.80;

            // Heston and Elices starting values (vertices) in vector form.  Add random increment about each starting value
            double[,] startH = new double[5,6];
            double[,] startE = new double[4,5];
            for(int j=0;j<=5;j++)
            {
                startH[0,j] = kappaS + BA.RandomNum(-0.05,0.05)*kappaS;
                startH[1,j] = thetaS + BA.RandomNum(-0.05,0.05)*thetaS;
                startH[2,j] = sigmaS + BA.RandomNum(-0.05,0.05)*sigmaS;
                startH[3,j] = rhoS   + BA.RandomNum(-0.05,0.05)*rhoS;
                startH[4,j] = v0S    + BA.RandomNum(-0.05,0.05)*v0S;
            }
            for(int j=0;j<=4;j++)
            {
                startE[0,j] = startH[0,j];
                startE[1,j] = startH[1,j];
                startE[2,j] = startH[2,j];
                startE[3,j] = startH[3,j];
            }

            // Settings for the Nelder Mead algorithm
            NMSet nmsettings = new NMSet();
            nmsettings.MaxIters = 50;		    // Maximum number of iterations
            nmsettings.Tolerance = 1e-20;		// Tolerance on best and worst function values
            nmsettings.ofsettingsE = ofsetE;    // Settings for the Elices objective function
            nmsettings.ofsettingsH = ofsetH;    // Settings for the Heston objective function

            // Heston Parameter Estimation all maturities ===========================================================================================
            nmsettings.choice = "Heston";
            int N = 5;
            nmsettings.N = N;
            double[] ParamH = new double[7];
            ParamH = NM.NelderMead(OF.f,nmsettings,startH);

            // Elices Parameter Estimation maturity by maturity ====================================================================================
            nmsettings.choice = "Elices";
            N = 4;
            nmsettings.N = N;
            double v0 = ParamH[4];
            ofsetE.v0 = v0;
            double[][] ParamTD = new double[4][];
            double[] B = new double[6];
            double[] LossFunction = new double[NT];

            // Initialize maturities
            List<double> MatsList = new List<double>();
            double[] Mats;
            MatsList.Add(tau[0]);
            MatsList.Add(tau[1]);
            Mats = MatsList.ToArray();

            // Loop through the maturities, estimating parameters
            for(int t=0;t<=NT-1;t++)
            {
                ofsetE.t = t;
                ofsetE.data.MktIV    = MktIV[t];            // Vector of market Implied vol
                ofsetE.data.MktPrice = MktPrice[t];         // Vector of market prices
                ofsetE.paramfixed = ParamTD;                // Fixed parameters
                if(t>=1)
                    for(int j=0;j<=N;j++)                       // Starting values are last period's estimates
                    {
                        startE[0,j] = B[0] + BA.RandomNum(-0.01,0.01)*B[0];
                        startE[1,j] = B[1] + BA.RandomNum(-0.01,0.01)*B[1];
                        startE[2,j] = B[2] + BA.RandomNum(-0.01,0.01)*B[2];
                        startE[3,j] = B[3] + BA.RandomNum(-0.01,0.01)*B[3];
                    }
                if(t>=2)
                {
                    MatsList.Add(tau[t]);
                    Mats = MatsList.ToArray();
                }
                ofsetE.T = Mats;
                nmsettings.ofsettingsE = ofsetE;
                B = NM.NelderMead(OF.f,nmsettings,startE);
                ParamTD[t] = new double[4] { B[0],B[1],B[2],B[3] };
                LossFunction[t] = B[4];
                Console.WriteLine("------ Finished maturity {0:F0} ------",ofsetE.t);
            }
            // Write the parameters
            Console.WriteLine("Heston statict parameters");
            Console.WriteLine("------------------------------------------------------------------");
            Console.WriteLine("    kappa      theta      sigma      rho         v0    LossFunction");
            Console.WriteLine("------------------------------------------------------------------");
            Console.WriteLine("{0,10:F5} {1,10:F5} {2,10:F5} {3,10:F5} {4,10:F5} {5,10:F5}",
                   ParamH[0],ParamH[1],ParamH[2],ParamH[3],ParamH[4],ParamH[5]);
            Console.WriteLine("------------------------------------------------------------------");
            Console.WriteLine(" ");
            Console.WriteLine("Elices Time dependent parameters");
            Console.WriteLine("-------------------------------------------------------------");
            Console.WriteLine("    kappa      theta      sigma      rho     LossFunction");
            Console.WriteLine("-------------------------------------------------------------");
            for(int t=0;t<=NT-1;t++)
                Console.WriteLine("{0,10:F5} {1,10:F5} {2,10:F5} {3,10:F5} {4,10:F5}",
                    ParamTD[t][0],ParamTD[t][1],ParamTD[t][2],ParamTD[t][3],LossFunction[t]);
            Console.WriteLine("-------------------------------------------------------------");


            // Fitting Elices prices and implied vols ==========================================================================================================
            // Matrices for prices, parameters, and fixed parameters
            double[,] PriceE = new double[NK,NT];
            double[] param = new double[4];
            double[][] paramfixed = new double[NT-1][];
            paramfixed[0] = ParamTD[0];

            // Define the dynamic array for the maturities and assign the first two maturities
            List<double> MatList = new List<double>();
            double[] Mat;
            MatList.Add(tau[0]);
            MatList.Add(tau[1]);
            Mat = MatList.ToArray();

            // Elices Implied volatilities
            double[,] IVE = new double[NK,NT];
            double a = 0.01;
            double b = 3.0;
            double tol = 1.0e-5;
            int MaxIter = 1000;
            double ElicesIVRMSE = 0.0;

            // Find the prices and implied vols
            for(int t=0;t<=NT-1;t++)
            {
                for(int j=0;j<=3;j++)
                    param[j] = ParamTD[t][j];
                if(t==0)
                    for(int k=0;k<=NK-1;k++)
                        PriceE[k,0] = EA.ElicesPrice(PutCall,opsettings.S,K[k],Mat,opsettings.r,opsettings.q,param,v0,opsettings.trap,X,W);
                if(t==1)
                    for(int k=0;k<=NK-1;k++)
                        PriceE[k,1] = EA.ElicesPrice(PutCall,opsettings.S,K[k],Mat,opsettings.r,opsettings.q,param,v0,opsettings.trap,X,W,paramfixed);
                if(t>=2)
                {
                    MatList.Add(tau[t]);
                    Mat = MatList.ToArray();
                    paramfixed[t-1] = ParamTD[t-1];
                    for(int k=0;k<=NK-1;k++)
                        PriceE[k,t] = EA.ElicesPrice(PutCall,opsettings.S,K[k],Mat,opsettings.r,opsettings.q,param,v0,opsettings.trap,X,W,paramfixed);
                }
                for(int k=0;k<=NK-1;k++)
                {
                    IVE[k,t] = BA.BisecBSIV(PutCall,opsettings.S,K[k],opsettings.r,opsettings.q,T[t],a,b,PriceE[k,t],tol,MaxIter);
                    ElicesIVRMSE += Math.Pow(IVE[k,t] - MktIV[t][k],2.0);
                }
            }
            // Fitting Heston prices and implied vols ==========================================================================================================
            double[,] PriceH = new double[NK,NT];
            double[,] IVH = new double[NK,NT];
            HParam ParamTI = new HParam();
            ParamTI.kappa = ParamH[0];
            ParamTI.sigma = ParamH[1];
            ParamTI.theta = ParamH[2];
            ParamTI.rho   = ParamH[3];
            ParamTI.v0    = ParamH[4];
            double HestonIVRMSE = 0.0;
            
            for (int k=0; k<=NK-1; k++)
                for(int t=0;t<=NT-1;t++)
                {
                    PriceH[k,t] = HP.HestonPriceGaussLaguerre(ParamTI,ofsetH.opsettings,K[k],T[t],PutCall,X,W);
                    IVH[k,t] = BA.BisecBSIV(PutCall,opsettings.S,K[k],opsettings.r,opsettings.q,T[t],a,b,PriceH[k,t],tol,MaxIter);
                    HestonIVRMSE += Math.Pow(IVH[k,t] - MktIVH[t,k],2.0);
                }
            Console.WriteLine("Heston IVRMSE {0,2:E5}",HestonIVRMSE);
            Console.WriteLine("Elices IVRMSE {0,2:E5}",ElicesIVRMSE);
        }
    }
}
