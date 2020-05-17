using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace Double_Heston_Simulation
{
    partial class DHSim
    {
        static void Main(string[] args)
        {
            // Classes
            HestonPriceDH HPDH = new HestonPriceDH();
            EulerAlfonsiSimulation EA = new EulerAlfonsiSimulation();
            QESimulation QE = new QESimulation();
            TVSimulation TV = new TVSimulation();

            // Spot price, risk free rate, dividend yield
            double S0 = 61.90;
            double K  = 61.90;
            double Mat = 1.0;
            double rf = 0.03;
            double q  = 0;
            string PutCall = "C";

            // Double Heston parameter values
            DHParam param;
            param.v01 = 0.36;
            param.v02 = 0.49;
            param.sigma1 =  0.1;
            param.sigma2 =  0.2;
            param.kappa1 =  0.9;
            param.kappa2 =  1.2;
            param.rho1   = -0.5;
            param.rho2   = -0.5;
            param.theta1 =  0.1;
            param.theta2 =  0.15;

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
            // Settings for the option
            OpSet settings;
            settings.S = S0;
            settings.K = K;
            settings.T = Mat;
            settings.r = rf;
            settings.q = q;
            settings.PutCall = PutCall;
            settings.trap = 1;

            // Closed form price
            double TruePrice = HPDH.DoubleHestonPriceGaussLaguerre(param,settings,X,W);

            // Simulation settings
            int NS = 10000;
            int NT = 1000;

            // Settings for the clock;
            Stopwatch sw = new Stopwatch();

            // Select the desired simulation scheme
            double SimPrice = 0.0;
            string scheme = "ZhuTV";

            // Simulation prices
            sw.Reset();
            sw.Start();
            if (scheme == "Euler")
                SimPrice = EA.DHEulerAlfonsiSim(scheme,param,S0,K,Mat,rf,q,NT,NS,PutCall);
            else if (scheme == "Alfonsi")
                SimPrice = EA.DHEulerAlfonsiSim(scheme,param,S0,K,Mat,rf,q,NT,NS,PutCall);
            else if (scheme == "ZhuEuler")
                SimPrice = TV.DHTransVolSim(scheme,param,S0,K,Mat,rf,q,NT,NS,PutCall);
            else if (scheme == "ZhuTV")
                SimPrice = TV.DHTransVolSim(scheme,param,S0,K,Mat,rf,q,NT,NS,PutCall);
            else if (scheme == "QE")
                SimPrice = QE.DHQuadExpSim(param,S0,K,Mat,rf,q,NT,NS,PutCall);
            sw.Stop();
            TimeSpan ts = sw.Elapsed;
            double Error = TruePrice - SimPrice;
            double ErrorP = Error/TruePrice*100;


            // Output the results
            Console.WriteLine("Using {0:0} simulations and {1:0} time steps",NS,NT);
            Console.WriteLine("Double Heston Simulation scheme : {0}",scheme);
            Console.WriteLine("--------------------------------------------------------------");
            Console.WriteLine("Scheme              Price    $Error    %Error   Sec   mSec");
            Console.WriteLine("--------------------------------------------------------------");
            Console.WriteLine("Closed form         {0,5:F4}",TruePrice);
            Console.WriteLine("Simulation          {0,5:F4} {1,8:F4} {2,8:F4} {3,5:F0} {4,5:F0}",SimPrice,Error,ErrorP,ts.Seconds,ts.Milliseconds);
            Console.WriteLine("--------------------------------------------------------------");
        }
    }
}
