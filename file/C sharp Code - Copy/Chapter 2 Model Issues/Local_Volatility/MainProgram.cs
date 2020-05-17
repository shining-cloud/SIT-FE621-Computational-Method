using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Numerics;

namespace Local_Volatility
{
    class SPX_LV
    {
        static void Main(string[] args)
        {
            // Gauss Laguerre 32 abscissas and weights
            double[] x = new Double[32];
            double[] w = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    x[k] = double.Parse(bits[0]);
                    w[k] = double.Parse(bits[1]);
                }

            // Size of strikes (NK) and maturities (NT)
            int NK = 27;
            int NT = 4;

            // Market SPX implied vol
            double[,] MktIV = new Double[NK,NT];
            using(TextReader reader = File.OpenText("../../SPX_MktIV.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    MktIV[k,0] = double.Parse(bits[0]);
                    MktIV[k,1] = double.Parse(bits[1]);
                    MktIV[k,2] = double.Parse(bits[2]);
                    MktIV[k,3] = double.Parse(bits[3]);
                }

            // Heston SPX implied vol
            double[,] HestonIV = new Double[NK,NT];
            using(TextReader reader = File.OpenText("../../SPX_HestonIV.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    HestonIV[k,0] = double.Parse(bits[0]);
                    HestonIV[k,1] = double.Parse(bits[1]);
                    HestonIV[k,2] = double.Parse(bits[2]);
                    HestonIV[k,3] = double.Parse(bits[3]);
                }

            // SPX strikes
            double[] K = new Double[NK];
            using(TextReader reader = File.OpenText("../../SPX_K.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    K[k] = double.Parse(bits[0]);
                }

            // SPX maturities
            double[] T = new Double[4] {0.019178082191781, 0.041095890410959, 0.117808219178082, 0.194520547945205};

            // SPX spot price, risk free rate, and dividend yield
            double S = 1164.97;
            double rf = 0.0;
            double q = 0.0;

            // Input the pre-calculated Heston parameters
            double v0     =  0.12581983608881;
            double theta  =  0.00000000008998;
            double kappa  =  0.02271544261755;
            double sigma  =  0.77411910754829;
            double rho    = -0.95555303183327;
            double lambda = 0.0;

            // Calculate the local volatilities
            LocalVol LV = new LocalVol();
            int trap = 1;
            double dt = 0.0001;
            double dK = 0.01;
            double[,] LVFD = new Double[NK,NT];
            double[,] LVAN = new Double[NK,NT];
            double[,] LVAP = new Double[NK,NT];
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<=NT-1;t++)
                {
                    LVFD[k,t] = LV.HestonLVFD(S,K[k],T[t],rf,q,kappa,theta,sigma,v0,rho,lambda,trap,x,w,dt,dK);
                    LVAN[k,t] = LV.HestonLVAnalytic(S,K[k],T[t],rf,q,kappa,theta,sigma,lambda,v0,rho,x,w,trap);
                    LVAP[k,t] = LV.HestonLVApprox(S,K[k],T[t],kappa,theta,sigma,v0,rho);
                }

            // Write the Finite Difference LV to a text file
            using(var writer = new StreamWriter("../../SPX_LocalVol_Finite_Differences.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    for(int t=0;t<=NT-1;t++)
                    {
                        writer.Write(LVFD[k,t]);
                        writer.Write(' ');
                    }
                    writer.WriteLine();
                }
            // Write the Analytic LV to a text file
            using(var writer = new StreamWriter("../../SPX_LocalVol_Analytic.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    for(int t=0;t<=NT-1;t++)
                    {
                        writer.Write(LVAN[k,t]);
                        writer.Write(' ');
                    }
                    writer.WriteLine();
                }
            // Write the approximate LV to a text file
            using(var writer = new StreamWriter("../../SPX_LocalVol_Approximate.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    for(int t=0;t<=NT-1;t++)
                    {
                        writer.Write(LVAP[k,t]);
                        writer.Write(' ');
                    }
                    writer.WriteLine();
                }

            // Output the results to the console
            Console.WriteLine("First Maturity --------------------------");
            Console.WriteLine("Approximate  Analytic   FiniteDifference ");
            Console.WriteLine("-----------------------------------------");
            for(int k=0;k<=NK-1;k++)
                Console.WriteLine("{0:F4} {1,12:F4} {2,12:F4}",LVAP[k,0],LVAN[k,0],LVFD[k,0]);

            Console.WriteLine(" ");
            Console.WriteLine("Second Maturity -------------------------");
            Console.WriteLine("Approximate  Analytic   FiniteDifference ");
            Console.WriteLine("-----------------------------------------");
            for(int k=0;k<=NK-1;k++)
                Console.WriteLine("{0:F4} {1,12:F4} {2,12:F4}",LVAP[k,1],LVAN[k,1],LVFD[k,1]);

            Console.WriteLine(" ");
            Console.WriteLine("Third Maturity --------------------------");
            Console.WriteLine("Approximate  Analytic   FiniteDifference ");
            Console.WriteLine("-----------------------------------------");
            for(int k=0;k<=NK-1;k++)
                Console.WriteLine("{0:F4} {1,12:F4} {2,12:F4}",LVAP[k,2],LVAN[k,2],LVFD[k,2]);

            Console.WriteLine(" ");
            Console.WriteLine("Fourth Maturity -------------------------");
            Console.WriteLine("Approximate  Analytic   FiniteDifference ");
            Console.WriteLine("-----------------------------------------");
            for(int k=0;k<=NK-1;k++)
                Console.WriteLine("{0:F4} {1,12:F4} {2,12:F4}",LVAP[k,3],LVAN[k,3],LVFD[k,3]);
            Console.WriteLine("-----------------------------------------");
            Console.WriteLine(" ");
        }
    }
}

