using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Price_Gauss_Laguerre_Consolidated
{
    class HestonPriceCon
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
            HParam param = new HParam();
            param.rho = -0.9;					// Heston Parameter: Correlation
            param.kappa = 2;					// Heston Parameter: Mean reversion speed
            param.theta = 0.05;				    // Heston Parameter: Mean reversion level 
            param.lambda = 0.0;				    // Heston Parameter: Risk parameter 
            param.sigma = 0.3;				    // Heston Parameter: Volatility of Variance
            param.v0 = 0.035;					// Heston Parameter: Current Variance

            OpSet settings = new OpSet();
            settings.S = 100.0;				    // Spot Price
            settings.K = 100.0;				    // Strike Price
            settings.T = 0.25;			        // Maturity in Years
            settings.r = 0.05;					// Interest Rate
            settings.q = 0.02;                  // Dividend yield
            settings.PutCall = "P";             // "P"ut or "C"all
            settings.trap = 1;                  // 1="Little Trap" characteristic function

            // The Heston price
            HestonPriceConsolidated HPC = new HestonPriceConsolidated();
            double Price = HPC.HestonPriceConsol(param,settings,x,w);
            Console.WriteLine("Heston price using consolidated integral");
            Console.WriteLine("------------------------------------------------ ");
            Console.WriteLine("Option Flavor =  {0,0:F5}",settings.PutCall);
            Console.WriteLine("Strike Price  =  {0,0:0}",settings.S);
            Console.WriteLine("Maturity      =  {0,0:F2}",settings.T);
            Console.WriteLine("Price         =  {0,0:F4}",Price);
            Console.WriteLine("------------------------------------------------ ");
            Console.WriteLine(" ");
        }
    }
}

