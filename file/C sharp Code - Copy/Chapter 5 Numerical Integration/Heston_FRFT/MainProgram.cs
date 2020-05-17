using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;


namespace Heston_FRFT
{
    class MainProgram
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] x = new Double[32];
            double[] w = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    x[k] = double.Parse(bits[0]);
                    w[k] = double.Parse(bits[1]);
                }
            }
            // Heston parameters
            HParam param = new HParam();
            param.kappa = 0.2;
            param.theta = 0.05;
            param.sigma = 0.3;
            param.v0 = 0.05;
            param.rho = -0.7;
            param.lambda = 0.0;

            // Option settings
            OpSet settings = new OpSet();
            settings.S = 50.0;
            settings.T = 0.5;
            settings.r = 0.03;
            settings.q = 0.05;
            settings.PutCall = "C";
            settings.trap = 1;

            // FFT settings
            double M = Math.Pow(2,10);
            int N = Convert.ToInt16(M);
            double alpha = 1.5;
            double uplimit = 100.0;
            string rule;
            
            // Set the integration increment
            double eta = 1e-2;

            // Set the log-strike increment
            double lambdainc = 1e-2;

            // Alternative: Set the log-strike increment so that the strike range is Spot +/ 5
            //double lambdainc = 2.0/M*Math.Log(settings.S/(settings.S - 5.0));

            // Obtain the FFT strikes and calls using the Trapezoidal rule
            HFRFT FRFT = new HFRFT();
            double[,] PriceTrapz = new double[N,2];
            rule = "Trapezoidal";
            PriceTrapz = FRFT.HestonFRFT(N,uplimit,alpha,rule,lambdainc,eta,param,settings);

            // Obtain the FFT strikes and calls using Simpson's rule
            double[,] PriceSimps = new double[N,2];
            rule = "Simpsons";
            PriceSimps = FRFT.HestonFRFT(N,uplimit,alpha,rule,lambdainc,eta,param,settings);

            // Obtain the extact price and compare to FFT price
            HestonPrice HP = new HestonPrice();
            double[] PriceHeston = new double[N];
            double[] ErrorTrapz = new double[N];
            double[] ErrorSimps = new double[N];
            for(int j=0;j<=N-1;j++)
            {
                settings.K = PriceTrapz[j,0];
                PriceHeston[j] = HP.HestonPriceGaussLaguerre(param,settings,x,w);
                ErrorTrapz[j] = (PriceTrapz[j,1] - PriceHeston[j])/PriceHeston[j]*100.0;
                ErrorSimps[j] = (PriceSimps[j,1] - PriceHeston[j])/PriceHeston[j]*100.0;
            }

            // Output the results near the money
            int start = N/2 - 3;
            int end   = N/2 + 3;
            Console.WriteLine("FRFT strike range from {0:F2} to {1:F2}",PriceTrapz[0,0],PriceTrapz[N-1,0]);
            Console.WriteLine();
            Console.WriteLine("Sample Strikes");
            Console.WriteLine("---------------------------------------------------------------------------");
            Console.WriteLine("Strike        Exact     FRFT Trapz  FRFT Simpson  %ErrorTrapz  %ErrorSimpson");
            Console.WriteLine("---------------------------------------------------------------------------");
            for(int j=start;j<=end;j++)
                Console.WriteLine("{0:F4} {1,12:F4} {2,12:F4} {3,12:F4} {4,12:F4} {5,12:F4}",PriceTrapz[j,0],PriceHeston[j],PriceTrapz[j,1],PriceSimps[j,1],ErrorTrapz[j],ErrorSimps[j]);
            Console.WriteLine("---------------------------------------------------------------------------");
            Console.WriteLine("Integration increment {0:F10}",eta);
            Console.WriteLine("Log strike increment {0:F10}",lambdainc);
        }
    }
}


