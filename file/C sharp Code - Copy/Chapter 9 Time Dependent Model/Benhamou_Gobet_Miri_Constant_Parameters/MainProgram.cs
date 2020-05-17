using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Benhamou_Gobet_Miri_Constant_Parameters
{
    public partial class BGM_Constant
    {
        static void Main(string[] args)
        {
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

            // Reproduces Table 3 and 4 of Benhamou, Gobet, and Miri
            //"Time Dependent Heston Model"
            // SIAM Journal on Financial Mathematics, Vol. 1, (2010)

            // Spot price, risk free rate, dividend yield
            OpSet settings = new OpSet();
            settings.S  = 100;
            settings.r = 0.0;
            settings.q  = 0.0;
            settings.PutCall = "C";
            settings.trap   = 1;       // Little trap formulation

            // Heston model parameters.
            // BGM Tables 3 and 4
            HParam param = new HParam();
            param.kappa  = 3.0;     // Volatility reversion speed
            param.theta  = 0.06;    // Volatility reversion level
            param.sigma  = 0.3;     // Volatility of variance
            param.rho    = 0.0;     // Correlation
            param.v0     = 0.04;    // Initial variance

            // Input the maturities and the strike prices
            double[,] T = new double[,] {
                {0.25,	0.25,	0.25,	0.25,	0.25,	0.25,	0.25,	0.25},
                {0.50,	0.50,	0.50,	0.50,	0.50,	0.50,	0.50,	0.50},
                {1.00,	1.00,	1.00,	1.00,	1.00,	1.00,	1.00,	1.00},
                {2.00,	2.00,	2.00,	2.00,	2.00,	2.00,	2.00,	2.00},
                {3.00,	3.00,	3.00,	3.00,	3.00,	3.00,	3.00,	3.00},
                {5.00,	5.00,	5.00,	5.00,	5.00,	5.00,	5.00,	5.00},
                {7.00,	7.00,	7.00,	7.00,	7.00,	7.00,	7.00,	7.00},
                {10.00,	10.00,	10.00,	10.00,	10.00,	10.00,	10.00,	10.00}};
            double[,] K = new double[,] {
                 {70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 125.0, 130},
	             {60.0, 70.0, 80.0, 100.0, 110.0, 130.0, 140.0, 150},
	             {50.0, 60.0, 80.0, 100.0, 120.0, 150.0, 170.0, 180},
	             {40.0, 50.0, 70.0, 100.0, 130.0, 180.0, 210.0, 240},
	             {30.0, 40.0, 60.0, 100.0, 140.0, 200.0, 250.0, 290},
	             {20.0, 30.0, 60.0, 100.0, 150.0, 250.0, 320.0, 400},
	             {10.0, 30.0, 50.0, 100.0, 170.0, 300.0, 410.0, 520},
	             {10.0, 20.0, 50.0, 100.0, 190.0, 370.0, 550.0, 730}};

            // Settings for the Bisection algorithm for implied vol
            double a = 1e-2;
            double b = 5.0;
            int MaxIter = 3000;
            double Tol = 1e-8;

            // Obtain the approximate and exact put prices, and their implied vols
            HestonPrice HP = new HestonPrice();
            BGMPrice BGM = new BGMPrice();
            BisectionAlgo BA = new BisectionAlgo();
            double[,] ApproxCall = new double[8,8];
            double[,] ExactCall  = new double[8,8];
            double[,] ApproxIV  = new double[8,8];
            double[,] ExactIV   = new double[8,8];
            for(int i=0;i<=7;i++)
                for(int j=0;j<=7;j++)
                {
                    ApproxCall[i,j] = BGM.BGMApproxPrice(param,settings,K[i,j],T[i,j]);
                    ExactCall[i,j]  = HP.HestonPriceGaussLaguerre(param,settings,K[i,j],T[i,j],X,W);
                    ApproxIV[i,j] = BA.BisecBSIV(settings,K[i,j],T[i,j],a,b,ApproxCall[i,j],Tol,MaxIter) * 100.0;
                    ExactIV[i,j]  = BA.BisecBSIV(settings,K[i,j],T[i,j],a,b,ExactCall[i,j],Tol,MaxIter)  * 100.0;
                }

            // Output Table 4
            Console.WriteLine("Benhamou, Gobet, Miri  Table 4");
            Console.WriteLine("Exact Call / Approximate Call");
            Console.WriteLine("------------------------------------------------------------------------");
            for(int i=0;i<=7;i++)
            {
                Console.WriteLine("{0,8:F2} {1,8:F2} {2,8:F2} {3,8:F2} {4,8:F2} {5,8:F2} {6,8:F2} {7,8:F2}",
                                   ExactCall[i,0],ExactCall[i,1],ExactCall[i,2],ExactCall[i,3],ExactCall[i,4],ExactCall[i,5],ExactCall[i,6],ExactCall[i,7]);
                Console.WriteLine("{0,8:F2} {1,8:F2} {2,8:F2} {3,8:F2} {4,8:F2} {5,8:F2} {6,8:F2} {7,8:F2}",
                                   ApproxCall[i,0],ApproxCall[i,1],ApproxCall[i,2],ApproxCall[i,3],ApproxCall[i,4],ApproxCall[i,5],ApproxCall[i,6],ApproxCall[i,7]);
            }
            Console.WriteLine("------------------------------------------------------------------------");

            // Output Table 3
            Console.WriteLine("Benhamou, Gobet, Miri  Table 3");
            Console.WriteLine("Exact Call IV / Approx Call IV");
            Console.WriteLine("------------------------------------------------------------------------");
            for(int i=0;i<=7;i++)
            {
                Console.WriteLine("{0,8:F2} {1,8:F2} {2,8:F2} {3,8:F2} {4,8:F2} {5,8:F2} {6,8:F2} {7,8:F2}",
                                   ExactIV[i,0],ExactIV[i,1],ExactIV[i,2],ExactIV[i,3],ExactIV[i,4],ExactIV[i,5],ExactIV[i,6],ExactIV[i,7]);
                Console.WriteLine("{0,8:F2} {1,8:F2} {2,8:F2} {3,8:F2} {4,8:F2} {5,8:F2} {6,8:F2} {7,8:F2}",
                                   ApproxIV[i,0],ApproxIV[i,1],ApproxIV[i,2],ApproxIV[i,3],ApproxIV[i,4],ApproxIV[i,5],ApproxIV[i,6],ApproxIV[i,7]);
            }
            Console.WriteLine("------------------------------------------------------------------------");

            // Output ATM values only
            Console.WriteLine("ATM Values only");
            Console.WriteLine("Exact IV       Approx IV   Exact Call  Approx Call");
            Console.WriteLine("----------------------------------------------------");
            for(int i=0;i<=7;i++)
                Console.WriteLine("{0,8:F2} {1,12:F2} {2,12:F2} {3,12:F2}",ExactIV[i,3],ApproxIV[i,3],ExactCall[i,3],ApproxCall[i,3]);
            Console.WriteLine("----------------------------------------------------");

        }
    }
}
