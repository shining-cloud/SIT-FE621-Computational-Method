using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Benhamou_Gobet_Miri_Piecewise_DIA_Estimation
{
    class BGM_Estimation
    {
        static void Main(string[] args)
        {
            // Classes
            BGMPrice BGM = new BGMPrice();
            BisectionAlgo BA = new BisectionAlgo();

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

            // DJIA ETF (ticker DIA) put data
            double[,] PutIV = new double[13,4] {
                {0.1962, 0.1947, 0.2019, 0.2115}, {0.1910, 0.1905, 0.1980, 0.2082},
                {0.1860, 0.1861, 0.1943, 0.2057}, {0.1810, 0.1812, 0.1907, 0.2021},
                {0.1761, 0.1774, 0.1871, 0.2000}, {0.1718, 0.1743, 0.1842, 0.1974},
                {0.1671, 0.1706, 0.1813, 0.1950}, {0.1644, 0.1671, 0.1783, 0.1927},
                {0.1645, 0.1641, 0.1760, 0.1899}, {0.1661, 0.1625, 0.1743, 0.1884},
                {0.1701, 0.1602, 0.1726, 0.1862}, {0.1755, 0.1610, 0.1716, 0.1846},
                {0.1796, 0.1657, 0.1724, 0.1842}};
            double[] K = new double[] { 124.0,125.0,126.0,127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0,136.0 };
            double[] T = new double[] { 37.0/365.0,72.0/365.0,135.0/365.0,226.0/365.0, };
            double[,] MktIV = PutIV;
            double S = 129.14;
            double r = 0.0010;
            double q  = 0.0068;
            int NT = 4;
            int NK = 13;

            // Create the market prices
            BlackScholesPrice BS = new BlackScholesPrice();
            string[,] PutCall = new string[NK,NT];
            double[,] MktPrice = new double[NK,NT];
            for(int t=0;t<NT;t++)
                for(int k=0;k<NK;k++)
                {
                    PutCall[k,t] = "P";
                    MktPrice[k,t] = BS.BlackScholes(S,K[k],T[t],r,q,MktIV[k,t],PutCall[k,t]);
                }

            // Bounds on the parameter estimates
            double e = 1e-3;
            double kappaU = 20.0;       double kappaL = e;
            double thetaU =  3.0;       double thetaL = e;
            double sigmaU = 10.0;       double sigmaL = e;
            double v0U    =  3.0;       double v0L    = e;
            double rhoU   = 0.99;       double rhoL   = -0.99;
            int N = 2 + 3*NT;
            double[] ub = new double[N];
            double[] lb = new double[N];
            for(int j=0;j<=N;j++)
            {
                ub[0] = kappaU; lb[0] = kappaL;
                ub[1] = v0U; lb[1] = v0L;
                for(int k=1;k<=NT;k++)
                {
                    ub[3*k-1] = thetaU; lb[3*k-1] = thetaL;
                    ub[3*k]   = sigmaU; lb[3*k-1] = sigmaL;
                    ub[3*k+1] = rhoU; lb[3*k+1] = rhoL;
                }
            }

            // Starting values for the parameters
            // Parameter order is [kappa v0 | theta[0] sigma[0] rho[0] | theta[1] sigma[1] rho[1] | etc...];
            double kappaS =  2.0;
            double thetaS =  0.1;
            double sigmaS =  1.2;
            double v0S    =  0.05;
            double rhoS   = -0.5;
            double[,] s = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                s[0,j] = kappaS + BS.RandomNum(-0.01,0.01)*kappaS;
                s[1,j] = v0S    + BS.RandomNum(-0.01,0.01)*v0S;
                for(int k=1;k<=NT;k++)
                {
                    s[3*k-1,j] = thetaS + BS.RandomNum(-0.01,0.01)*thetaS;
                    s[3*k,j]   = sigmaS + BS.RandomNum(-0.01,0.01)*sigmaS;
                    s[3*k+1,j] = rhoS   + BS.RandomNum(-0.01,0.01)*rhoS;
                }
            }

            // Structure for the option settings
            OPSet opset;
            opset.S = S;
            opset.r = r;
            opset.q  = q;
            opset.trap = 1;
            opset.PutCall = "P";

            // Structure for the market data
            MktData data;
            data.MktIV = MktIV;
            data.MktPrice = MktPrice;
            data.K = K;
            data.T = T;
            data.PutCall = PutCall;

            // Structure for the objective function settings
            OFSet ofset;
            ofset.opset = opset;
            ofset.data = data;
            ofset.X = X;
            ofset.W = W;
            ofset.ObjFunction = 4;
            ofset.lb = lb;
            ofset.ub = ub;

            // Structure for the Nelder Mead algorithm settings
            NMSet nmset;
            nmset.ofset = ofset;
            nmset.MaxIters = 5000;
            nmset.Tolerance = 1e-8;
            nmset.N = N;

            // Run the Nelder Mead Algorithm
            NelderMeadAlgo NM = new NelderMeadAlgo();
            ObjectiveFunction OF = new ObjectiveFunction();
            double[] B = NM.NelderMead(OF.f,nmset,s);
            int NB = B.Length;

            // Separate the parameter vector into different vectors
            double kappa = B[0];
            double v0 = B[1];
            double[] theta = new double[NT];
            double[] sigma = new double[NT];
            double[] rho   = new double[NT];
            for(int k=0;k<=NT-1;k++)
            {
                theta[k] = B[3*k+2];
                sigma[k] = B[3*k+3];
                rho[k]   = B[3*k+4];
            }

            // Bisection algorithm settings
            double a = 0.01;
            double b = 4.0;
            double Tol = 1e-5;
            int MaxIter = 2500;
            
            // Find the implied volatilities
            double[,] ModelPrice = new double[NK,NT];
            double[,] ModelIV    = new double[NK,NT];
            double Error = 0.00;

            List<double> MatList   = new List<double>();
            List<double> thetaList = new List<double>();
            List<double> sigmaList = new List<double>();
            List<double> rhoList   = new List<double>();

            double[] Mat,Theta,Sigma,Rho;
            for(int t=0;t<NT;t++)
            {
                // Stack the parameters
                MatList.Add(T[t]);
                thetaList.Add(theta[t]);
                sigmaList.Add(sigma[t]);
                rhoList.Add(rho[t]);
                // Convert to arrays
                Mat = MatList.ToArray();
                Theta = thetaList.ToArray();
                Sigma = sigmaList.ToArray();
                Rho   = rhoList.ToArray();
                for(int k=0;k<NK;k++)
                {
                    ModelPrice[k,t] = BGM.BGMApproxPriceTD(kappa,v0,Theta,Sigma,Rho,opset,K[k],Mat);
                    ModelIV[k,t] = BA.BisecBSIV(opset,K[k],T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                    Error += Math.Pow(ModelIV[k,t] - MktIV[k,t],2.0) / Convert.ToDouble(NT*NK);
                }
            }

            // Output the results
            Console.WriteLine("----------------------------------------------------------");
            Console.WriteLine("    kappa      v0         theta      sigma        rho");
            Console.WriteLine("----------------------------------------------------------");
            for(int k=0;k<=NT-1;k++)
                Console.WriteLine("{0,10:F5} {1,10:F5} {2,10:F5} {3,10:F5} {4,10:F5}",kappa,v0,theta[k],sigma[k],rho[k]);

            Console.WriteLine("----------------------------------------------------------");
            Console.WriteLine("Value of objective function   {0:F5} ",B[NB-2]);
            Console.WriteLine("Number of iterations required {0:0}  ",B[NB-1]);
            Console.WriteLine("----------------------------------------------------------");
            Console.WriteLine(" ");
            Console.WriteLine("Model/Market implied volatilities");
            Console.WriteLine("----------------------------------------------------------");
            Console.WriteLine("Strike    Mat1       Mat2        Mat3       Mat4");
            Console.WriteLine("----------------------------------------------------------");
            for(int k=0;k<=NK-1;k++)
            {
                Console.WriteLine("{0,5:0} {1,10:F5} {2,10:F5} {3,10:F5} {4,10:F5}",K[k],ModelIV[k,0],ModelIV[k,1],ModelIV[k,2],ModelIV[k,3]);
                Console.WriteLine("{0,5:0} {1,10:F5} {2,10:F5} {3,10:F5} {4,10:F5}",K[k],MktIV[k,0],MktIV[k,1],MktIV[k,2],MktIV[k,3]);
                Console.WriteLine("----------------------------------------------------------");

            }
            Console.WriteLine(" ");
            Console.WriteLine("IV Estimation error {0,5:E5}", Error);
            Console.WriteLine(" ");
        }
    }
}

