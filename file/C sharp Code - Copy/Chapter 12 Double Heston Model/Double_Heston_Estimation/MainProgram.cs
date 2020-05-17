using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;
using System.Diagnostics;

namespace Double_Heston_Estimation
{
    public partial class DoubleHestonEstimation
    {
        static void Main(string[] args)
        {
            BlackScholesPrice BA = new BlackScholesPrice();
            NelderMeadAlgo NM = new NelderMeadAlgo();
            ObjectiveFunction OF = new ObjectiveFunction();
            ObjectiveFunctionSVC OFSVC = new ObjectiveFunctionSVC();

            // Settings for the Nelder Mead algorithm
            int N = 10;						// Number of Heston parameters 
            int NumIters = 1;			    // First Iteration
            int MaxIters = 5000;	        // Maximum number of iterations
            double Tolerance = 1e-3;		// Tolerance on best and worst function values
            int LossFunction = 3;           // Choice of loss function
                                            // 1=MSE, 2=RMSE, 3=IVMSE, 4=Christoffersen et al.

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
            int NK = 13;
            int NT = 4;
            double[,] MktIV = new double[13,4] {
                {0.1962, 0.1947, 0.2019, 0.2115}, {0.1910, 0.1905, 0.1980, 0.2082},
                {0.1860, 0.1861, 0.1943, 0.2057}, {0.1810, 0.1812, 0.1907, 0.2021},
                {0.1761, 0.1764, 0.1871, 0.2000}, {0.1718, 0.1743, 0.1842, 0.1974},
                {0.1671, 0.1706, 0.1813, 0.1950}, {0.1644, 0.1671, 0.1783, 0.1927},
                {0.1661, 0.1641, 0.1760, 0.1899}, {0.1661, 0.1625, 0.1743, 0.1884},
                {0.1701, 0.1602, 0.1726, 0.1871}, {0.1755, 0.1610, 0.1716, 0.1846},
                {0.1786, 0.1657, 0.1724, 0.1842}};
            double[] K = new double[13] { 124.0,125.0,126.0,127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0,136.0 };
            double[] T = new double[4] { 37.0/365.0,72.0/365.0,135.0/365.0,226.0/365.0 };
            string[,] PutCall = new String[NK,NT];
            for(int k=0;k<=NK-1;k++)
            {
                for(int t=0;t<=NT-1;t++)
                    PutCall[k,t] = "C";
            }

            // Option settings
            OpSet settings = new OpSet();
            settings.S = 129.14;
            settings.r = 0.0010;
            settings.q = 0.0068;
            settings.trap = 1;

            settings.r = 0.0;
            settings.q = 0.0;

            // Obtain the market prices
            double[,] MktPrice = new Double[NK,NT];
            for(int k=0;k<=NK-1;k++)
            {
                for(int t=0;t<=NT-1;t++)
                    MktPrice[k,t] = BA.BlackScholes(settings.S,K[k],T[t],settings.r,settings.q,MktIV[k,t],PutCall[k,t]);
            }

            // Market data
            MktData data = new MktData();
            data.MktIV = MktIV;
            data.MktPrice = MktPrice;
            data.K = K;
            data.T = T;
            data.PutCall = PutCall;

            // starting values for the Nelder Mean algorithm
            double kappaS =  3.0;
            double thetaS =  0.06;
            double sigmaS =  0.5;
            double v0S    =  0.04;
            double rhoS   = -0.8;

            // Create the vertices for the starting values
            double[,] s = new double[N,N+1];
            for(int j=0;j<=N;j++)
            {
                s[0,j] = kappaS + BA.RandomNum(-0.10,0.10);
                s[1,j] = thetaS + BA.RandomNum(-0.01,0.01);
                s[2,j] = sigmaS + BA.RandomNum(-0.01,0.01);
                s[3,j] = v0S    + BA.RandomNum(-0.01,0.01);
                s[4,j] = rhoS   + BA.RandomNum(-0.05,0.05);
                s[5,j] = kappaS + BA.RandomNum(-0.10,0.10);
                s[6,j] = thetaS + BA.RandomNum(-0.01,0.01);
                s[7,j] = sigmaS + BA.RandomNum(-0.01,0.01);
                s[8,j] = v0S    + BA.RandomNum(-0.01,0.01);
                s[9,j] = rhoS   + BA.RandomNum(-0.05,0.05);
            }

            // Upper and lower parameter bounds
            // kappa theta sigma v0 rho
            double e = 1e-5;
            double[] lb = new double[5] { e,e,e,e,-0.999 };
            double[] ub = new double[5] { 20.0,2.0,2.0,3.0,0.999 };

            // Settings for the clock;
            Stopwatch sw = new Stopwatch();

            // Estimation of parameters by Nelder Mean -- Ordinary method
            sw.Reset();
            sw.Start();
            double[] B = NM.NelderMead(OF.DHObjFun,N,NumIters,MaxIters,Tolerance,s,settings,data,X,W,LossFunction,lb,ub);
            sw.Stop();
            TimeSpan tsOrd = sw.Elapsed;

            // Estimation of parameters by Nelder Mean -- Strike Vector Computation method
            sw.Reset();
            sw.Start();
            double[] C = NM.NelderMead(OFSVC.DHObjFunSVC,N,NumIters,MaxIters,Tolerance,s,settings,data,X,W,LossFunction,lb,ub);
            sw.Stop();
            TimeSpan tsSVC = sw.Elapsed;

            // Output the results
            Console.WriteLine("---------------------------------------------------");
            Console.WriteLine("Double Heston parameter estimation");
            Console.WriteLine("Parameter            Ordinary method    SVC Method ");
            Console.WriteLine("---------------------------------------------------");
            Console.WriteLine("kappa1     {0,20:F5}      {1,10:F5} ",B[0],C[0]);
            Console.WriteLine("theta1     {0,20:F5}      {1,10:F5} ",B[1],C[1]);
            Console.WriteLine("sigma1     {0,20:F5}      {1,10:F5} ",B[2],C[2]);
            Console.WriteLine("v01        {0,20:F5}      {1,10:F5} ",B[3],C[3]);
            Console.WriteLine("rho1       {0,20:F5}      {1,10:F5} ",B[4],C[4]);
            Console.WriteLine(" ");
            Console.WriteLine("kappa2     {0,20:F5}      {1,10:F5} ",B[5],C[5]);
            Console.WriteLine("theta2     {0,20:F5}      {1,10:F5} ",B[6],C[6]);
            Console.WriteLine("sigma2     {0,20:F5}      {1,10:F5} ",B[7],C[7]);
            Console.WriteLine("v02        {0,20:F5}      {1,10:F5} ",B[8],C[8]);
            Console.WriteLine("rho2       {0,20:F5}      {1,10:F5} ",B[9],C[9]);
            Console.WriteLine("---------------------------------------------------");
            Console.WriteLine("Value of obj function  {0,7:F2}      {1,10:F2} ",B[10],C[10]);
            Console.WriteLine("Number of iterations   {0,6:0}        {1,6:0}  ",B[11],C[11]);
            Console.WriteLine("---------------------------------------------------");
            Console.WriteLine(" ");
            Console.WriteLine("Estimation time          Min   Sec  mSec");
            Console.WriteLine("-------------------------------------------");
            Console.WriteLine("Ordinary ObjFun {0,10:F0} {1,6:F0} {2,6:F0}",tsOrd.Minutes,tsOrd.Seconds,tsOrd.Milliseconds);
            Console.WriteLine("SVC ObjFun      {0,10:F0} {1,6:F0} {2,6:F0}",tsSVC.Minutes,tsSVC.Seconds,tsSVC.Milliseconds);
            Console.WriteLine("-------------------------------------------");
            Console.WriteLine(" ");
        }
    }
}
