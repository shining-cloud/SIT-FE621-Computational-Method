using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Mikhailov_and_Nogel_Estimation_DJIA
{
    class MikhailovNogelEstimation
    {
        static void Main(string[] args)
        {
            // Classes
            BlackScholesPrice BS = new BlackScholesPrice();
            BisectionAlgo BA = new BisectionAlgo();
            NelderMeadAlgo NM = new NelderMeadAlgo();
            ObjectiveFunction OF = new ObjectiveFunction();
            HestonPriceTD HPTD = new HestonPriceTD();

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
            int NT = 4;
            int NK = 13;
            double[,] PutIV = new double[13,4]{
                {0.1962, 0.1947, 0.2019, 0.2115}, {0.1910, 0.1905, 0.1980, 0.2082},
                {0.1860, 0.1861, 0.1943, 0.2057}, {0.1810, 0.1812, 0.1907, 0.2021},
                {0.1761, 0.1774, 0.1871, 0.2000}, {0.1718, 0.1743, 0.1842, 0.1974},
                {0.1671, 0.1706, 0.1813, 0.1950}, {0.1644, 0.1671, 0.1783, 0.1927},
                {0.1645, 0.1641, 0.1760, 0.1899}, {0.1661, 0.1625, 0.1743, 0.1884},
                {0.1701, 0.1602, 0.1726, 0.1862}, {0.1755, 0.1610, 0.1716, 0.1846},
                {0.1796, 0.1657, 0.1724, 0.1842}};
            double[] K = new double[13]{124.0, 125.0, 126.0, 127.0, 128.0, 129.0, 130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0,};
            double[] T = new double[4]{0.1014, 0.1973, 0.3699, 0.6192};
            string[,] PutCall = new string[13,4];
            double[,] MktPrice = new double[13,4];
            double Spot = 129.14;
            double r = 0.0010;
            double q = 0.0068;
            int trap = 1;
            int LossFunction = 2;

            for(int t=0;t<=NT-1;t++)
                for(int k=0;k<=NK-1;k++)
                {
                    PutCall[k,t] = "P";
                    MktPrice[k,t] = BS.BlackScholes(Spot,K[k],T[t],r,q,PutIV[k,t],PutCall[k,t]);
                }

            // Create the maturity increments
            double[] tau = new double[4];
            tau[0] = T[0];
            for (int t=1; t<=NT-1; t++)
	            tau[t] = T[t] - T[t-1];

            // Starting values and upper bounds
            double e = 1e-5;
            double[] lb = new double[5] {e,e,e,e,-0.999};               // Lower bound on the estimates
            double[] ub = new double[5] {20.0,2.0,2.0,3.0,0.999};       // Upper bound on the estimates
            
            // Starting values (vertices) in vector form.  Add random increment about each starting value
            int N = 5;
            double kappaS =  4.0;
            double thetaS =  0.1;
            double sigmaS =  1.5;
            double v0S    =  0.04;
            double rhoS   = -0.30;
            double[,] s = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                s[0,j] = kappaS + BA.RandomNum(-0.01,0.01);
                s[1,j] = thetaS + BA.RandomNum(-0.01,0.01);
                s[2,j] = sigmaS + BA.RandomNum(-0.01,0.01);
                s[3,j] = v0S    + BA.RandomNum(-0.01,0.01);
                s[4,j] = rhoS   + BA.RandomNum(-0.01,0.01);
            }

            // Arrays for old maturities and parameters, and current parameters (old parameters, but stacked)
            ArrayList OldTau = new ArrayList();
            double[,] param0 = new double[NT-1,5];
            double[,] paramTD = new double[NT,5];

            // Nelder Mead settings
            int MaxIters = 1000;
            double Tolerance = 1e-5;

            // Market data
            MktData data = new MktData();
            double[] MKIV = new double[NK];
            double[] MKPR = new double[NK];
            string[] PC = new string[NK];

            // Step-by-step parameter estimates
            double[] B = new double[6];
            double[] ObjectiveFun = new double[NT];
            double[] NumIteration = new double[NT];

            // First maturity parameter estimates ===================================================================
            double[] tau0 = {0.0};
            for(int k=0;k<=NK-1;k++)
            {
                MKIV[k] = PutIV[k,0];
                MKPR[k] = MktPrice[k,0];
                PC[k]  = PutCall[k,0];
            }
            data.MktIV = MKIV;
            data.PutCall = PC;
            data.MktPrice = MKPR;
            B = NM.NelderMead(OF.f,N,MaxIters,Tolerance,s,tau[0],tau0,param0,data,Spot,K,r,q,trap,X,W,LossFunction,lb,ub);
            Console.WriteLine("Finished parameter estimate set 1");
            for(int k=0;k<=4;k++)
                paramTD[0,k] = B[k];
            ObjectiveFun[0] = B[5];
            NumIteration[0] = B[6];

            // Remaining Maturity estimates and updated starting values ==================================================
            for(int mat=1;mat<=NT-1;mat++)
            {
                OldTau.Add(tau[mat-1]);
                for(int k=0;k<=4;k++)
                    param0[mat-1,k] = B[k];
                for(int k=0;k<=NK-1;k++)
                {
                    MKIV[k] = PutIV[k,mat];
                    MKPR[k] = MktPrice[k,mat];
                    PC[k]  = PutCall[k,mat];
                }
                data.MktIV = MKIV;
                data.PutCall = PC;
                data.MktPrice = MKPR;
                for(int j=0;j<=N;j++)
                {
                    s[0,j] = B[0] + BA.RandomNum(-0.01,0.01);
                    s[1,j] = B[1] + BA.RandomNum(-0.01,0.01);
                    s[2,j] = B[2] + BA.RandomNum(-0.01,0.01);
                    s[3,j] = B[3] + BA.RandomNum(-0.01,0.01);
                    s[4,j] = B[4] + BA.RandomNum(-0.01,0.01);
                }
                B = NM.NelderMead(OF.f,N,MaxIters,Tolerance,s,tau[mat],tau0,param0,data,Spot,K,r,q,trap,X,W,LossFunction,lb,ub);
                Console.WriteLine("Finished parameter estimate set {0}",mat+1);
                for(int k=0;k<=4;k++)
                    paramTD[mat,k] = B[k];
                ObjectiveFun[mat] = B[5];
                NumIteration[mat] = B[6];
            }

            // Fit the prices and implied volatilities ====================================================================
            // Bisection algorithm settings
            double a = 0.001;
            double b = 5.0;
            double Tol = 1e-5;
            int MaxIter = 5000;
            
            ArrayList OldTau00 = new ArrayList();
            double[] tau00 = { 0.0 };
            double[,] ModelPrice = new double[NK,NT];
            double[,] ModelIV = new double[NK,NT];
            Array.Clear(param0,0,15);

            // First maturity
            HParam param = new HParam();
            param.kappa = paramTD[0,0];
            param.theta = paramTD[0,1];
            param.sigma = paramTD[0,2];
            param.v0    = paramTD[0,3];
            param.rho   = paramTD[0,4];
            for(int k=0;k<=NK-1;k++)
            {
                ModelPrice[k,0] = HPTD.MNPriceGaussLaguerre(param,param0,tau[0],tau00,Spot,K[k],r,q,PutCall[k,0],trap,X,W);
                ModelIV[k,0] = BA.BisecBSIV(PutCall[k,0],Spot,K[k],r,q,T[0],a,b,ModelPrice[k,0],Tol,MaxIter);
            }

            // Remaining maturity
            for(int t=1;t<=NT-1;t++)
            {
                OldTau00.Add(tau[t-1]);
                param.kappa = paramTD[t,0];
                param.theta = paramTD[t,1];
                param.sigma = paramTD[t,2];
                param.v0    = paramTD[t,3];
                param.rho   = paramTD[t,4];
                for(int j=0;j<=4;j++)
                    param0[t-1,j] = paramTD[t-1,j];
                for(int k=0;k<=NK-1;k++)
                {
                    ModelPrice[k,t] = HPTD.MNPriceGaussLaguerre(param,param0,tau[t],tau00,Spot,K[k],r,q,PutCall[k,t],trap,X,W);
                    ModelIV[k,t] = BA.BisecBSIV(PutCall[k,t],Spot,K[k],r,q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                }
            }

            // Mean Square Implied Volatility Error
            double Error = 0.0;
            for (int t=0;t<=NT-1; t++)
                for (int k=0; k<=NK-1; k++)
                    Error += Math.Pow(ModelIV[k,t] - PutIV[k,t],2);

            Error /=  (Convert.ToDouble(NT)*Convert.ToDouble(NK));

            // Output the estimation result
            Console.WriteLine("-----------------------------------------------------------------------------");
            Console.WriteLine("Market/Model implied volatilities");
            Console.WriteLine("    Mat1       Mat2       Mat3       Mat4");
            for(int k=0;k<=NK-1;k++)
            {
                Console.WriteLine("--------------------------------------------- Strike {0} Market/Model",K[k]);
                Console.WriteLine("{0,10:F4} {1,10:F4} {2,10:F4} {3,10:F4}",PutIV[k,0],PutIV[k,1],PutIV[k,2],PutIV[k,3]);
                Console.WriteLine("{0,10:F4} {1,10:F4} {2,10:F4} {3,10:F4}",ModelIV[k,0],ModelIV[k,1],ModelIV[k,2],ModelIV[k,3]);
            }
            Console.WriteLine("-----------------------------------------------------------------------------");
            Console.WriteLine("Model Price ");
            for(int k=0;k<=NK-1;k++)
                Console.WriteLine("{0,10:0} {1,10:F4} {2,10:F4} {3,10:F4} {4,10:F4}",K[k],ModelPrice[k,0],ModelPrice[k,1],ModelPrice[k,2],ModelPrice[k,3]);
            Console.WriteLine("-----------------------------------------------------------------------------");
            Console.WriteLine("Results of time dependent parameter estimation");
            Console.WriteLine("  ");
            Console.WriteLine("Maturity   kappa   theta   sigma   v0       rho     ObjecFun  NumIterations");
            Console.WriteLine("-----------------------------------------------------------------------------");
            for(int t=0;t<=NT-1;t++)
                Console.WriteLine("{0} {1,10:F4} {2,7:F4} {3,7:F4} {4,7:F4} {5,8:F4} {6,8:F4} {7,7:0}",
                                  T[t],paramTD[t,0],paramTD[t,1],paramTD[t,2],paramTD[t,3],paramTD[t,4],ObjectiveFun[t],NumIteration[t]);
            Console.WriteLine("-----------------------------------------------------------------------------");
            Console.WriteLine("IVMSE from time dependent estimates {0,10:E5}", Error);
            Console.WriteLine("-----------------------------------------------------------------------------");
        }
    }
}

