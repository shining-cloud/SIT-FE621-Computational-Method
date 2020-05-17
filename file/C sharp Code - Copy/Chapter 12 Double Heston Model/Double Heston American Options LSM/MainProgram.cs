using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace Double_Heston_American_Options_LSM
{
    class DHAmerican
    {
        static void Main(string[] args)
        {
            // Classes
            HestonPriceDH HPDH = new HestonPriceDH();
            EulerAlfonsiSimulation EA = new EulerAlfonsiSimulation();
            TVSimulation TV = new TVSimulation();
            QESimulation QE = new QESimulation();
            Regression RE = new Regression();
            LSM LSM = new LSM();

            // Spot price, risk free rate, dividend yield
            double S0 = 61.90;
            double K  = 61.90;
            double Mat = 1.0;
            double rf = 0.03;
            double q  = 0.0;
            string PutCall = "P";

            // Double Heston parameter values
            DHParam param;
            param.v01 = 0.36;
            param.v02 = 0.49;
            param.sigma1 =  0.10;
            param.sigma2 =  0.20;
            param.kappa1 =  0.90;
            param.kappa2 =  1.20;
            param.rho1   = -0.5;
            param.rho2   = -0.5;
            param.theta1 =  0.10;
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
            double TrueEuroPrice = HPDH.DoubleHestonPriceGaussLaguerre(param,settings,X,W);

            // Simulation settings
            int NS = 50000;
            int NT = 1000;

            // Settings for the clock;
            Stopwatch sw = new Stopwatch();

            // Select the simulation scheme
            string scheme = "Euler";
            double LSMEuro,LSMAmer,Premium,CVPrice;

            // Zhu-Euler Double Heston prices
            DHSim Soutput = new DHSim();
            sw.Reset();
            sw.Start();
            if((scheme == "Euler") || (scheme == "Alfonsi"))
                Soutput = EA.DHEulerAlfonsiSim(scheme,param,S0,K,Mat,rf,q,NT,NS,PutCall);
            else if((scheme=="ZhuEuler") || (scheme=="ZhuTV"))
                Soutput = TV.DHTransVolSim(scheme,param,S0,K,Mat,rf,q,NT,NS,PutCall);
            else if(scheme == "QE")
                Soutput = QE.DHQuadExpSim(param,S0,K,Mat,rf,q,NT,NS,PutCall);
            double[,] S = Soutput.S;
            double SimEuro = Soutput.EuroPrice;
            sw.Stop();
            TimeSpan ts = sw.Elapsed;

            // LSM algorithm
            double[] LSMoutput = new double[2];
            double[,] St1 = RE.MTrans(S);
            LSMoutput = LSM.HestonLSM(St1,K,rf,q,Mat,NT,NS,PutCall);
            LSMEuro = LSMoutput[0];
            LSMAmer = LSMoutput[1];
            Premium = LSMAmer - LSMEuro;

            // Control variate American price
            CVPrice = TrueEuroPrice + Premium;

            // Output the results
            Console.WriteLine("----------------------------------------------------------------------");
            Console.WriteLine("Number of simulations {0:F0}",NS);
            Console.WriteLine("Number of time steps  {0:F0}",NT);
            Console.WriteLine("Simulation Scheme : {0}",scheme);
            Console.WriteLine("----------------------------------------------------------------------");
            Console.WriteLine("Method            Euro       Amer      Premium   ControlVar  Sec  mSec");
            Console.WriteLine("----------------------------------------------------------------------");
            Console.WriteLine("Closed        {0,10:F4}",TrueEuroPrice);
            Console.WriteLine("Simulation    {0,10:F4} {1,10:F4} {2,10:F4} {3,10:F4} {4,5:F0} {5,5:F0}",LSMEuro,LSMAmer,Premium,CVPrice,ts.Seconds,ts.Milliseconds);
            Console.WriteLine("----------------------------------------------------------------------");
        }
    }
}
