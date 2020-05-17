using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;
using System.Diagnostics;

namespace Heston_LSM
{
    public partial class LSM
    {
        static void Main(string[] args)
        {
            Stopwatch sw = new Stopwatch();
            TimeSpan ts = sw.Elapsed;
            
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

            // Clarke and Parrott settings
            double[] TruePrice = new double[5] { 2.0, 1.1070641,0.520030,0.213668,0.082036 };
            double[] Spot = new double[5] { 8.0, 9.0,10.0,11.0,12.0 };

            HParam param = new HParam();
            param.kappa = 5.0;
            param.theta = 0.16;
            param.sigma = 0.9;
            param.v0 = 0.0625;
            param.rho = 0.1;
            param.lambda = 0.0;

            OpSet settings = new OpSet();
            settings.K = 10.0;
            settings.T = 0.25;
            settings.r = 0.1;
            settings.q = 0.0;
            settings.PutCall = "P";
            settings.trap = 1;

            // Simulation settings
            int NT = 100;
            int NS = 5000;

            // Generate the correlated random variables
            RandomNumber RN = new RandomNumber();
            double[,] Zv = new double[NT,NS];
            double[,] Zs = new double[NT,NS];
            double rho = param.rho;
            sw.Reset();
            sw.Start();
            for(int t=0;t<=NT-1;t++)
                for(int s=0;s<=NS-1;s++)
                {
                    Zv[t,s] = RN.RandomNorm();
                    Zs[t,s] = rho*Zv[t,s] + Math.Sqrt(1-rho*rho)*RN.RandomNorm();
                }

            // Generate LSM prices using the correlated random variables generated above
            int M = Spot.Length;
            double[] ClosedEuro = new double[M];
            double[] LSMEuro = new double[M];
            double[] LSMAmer = new double[M];
            double[] CVAmer = new double[M];
            double[] ErrorLSM  = new double[M];
            double[] ErrorCV  = new double[M];
            double[] ErrorE  = new double[M];
            double TotalLSM=0.0,TotalCV=0.0,TotalE = 0.0;
            double[] LSMPrice = new double[2];

            // Obtain the prices
            LSMonteCarlo LS = new LSMonteCarlo();
            HestonPrice HP = new HestonPrice();
            Regression RE = new Regression();
            MMSimulation MM = new MMSimulation();
            for(int k=0;k<=M-1;k++)
            {
                // LSM Euro and American prices
                settings.S = Spot[k];
                LSMPrice = LS.HestonLSM(RE.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                LSMEuro[k] = LSMPrice[0];
                LSMAmer[k] = LSMPrice[1];
                // Closed Euro price
                ClosedEuro[k] = HP.HestonPriceGaussLaguerre(param,settings,X,W);
                // Control variate price
                CVAmer[k] = ClosedEuro[k] + (LSMAmer[k] - LSMEuro[k]);
                // The errors and total absolute errors
                ErrorLSM[k] = TruePrice[k] - LSMAmer[k];
                ErrorCV[k] = TruePrice[k] - CVAmer[k];
                ErrorE[k] = ClosedEuro[k] - LSMEuro[k];
                TotalLSM += Math.Abs(ErrorLSM[k]);
                TotalCV  += Math.Abs(ErrorCV[k]);
                TotalE   += Math.Abs(ErrorE[k]);
            }
            sw.Stop();
            ts = sw.Elapsed;
            
            Console.WriteLine("Comparison of Clarke and Parrott American Put prices with L-S Monte Carlo");
            Console.WriteLine("LSM uses {0:0} time steps and {1:0} stock paths",NT,NS);
            Console.WriteLine("-----------------------------------------------------------");
            Console.WriteLine(" S(0)  TrueAmer  LSMAmer    CVAmer    ClosedEuro  LSMEuro");
            Console.WriteLine("-----------------------------------------------------------");
            for(int k=0;k<=M-1;k++)
                Console.WriteLine("{0,3} {1,10:F6} {2,10:F6} {3,10:F6} {4,10:F6} {5,10:F6}",
                    Spot[k],TruePrice[k],LSMAmer[k],CVAmer[k],ClosedEuro[k],LSMEuro[k]);
            Console.WriteLine("-----------------------------------------------------------");
            Console.WriteLine("AbsError {0,15:F5} {1,10:F5} {2,21:F5}",TotalLSM,TotalCV,TotalE);
            Console.WriteLine("-----------------------------------------------------------");
            Console.WriteLine("Simulation time {0:0} minutes and {1,5:F3} seconds ",ts.Minutes,ts.Seconds);
        }
    }
}
