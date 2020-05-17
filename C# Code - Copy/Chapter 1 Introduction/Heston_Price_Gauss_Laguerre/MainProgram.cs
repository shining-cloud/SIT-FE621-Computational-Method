using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Price_Gauss_Laguerre
{
    class HestonPriceGaussLaguerre
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
            param.kappa = 1.5;
            param.theta = 0.04;
            param.sigma = 0.3;
            param.v0 = 0.05412;
            param.rho = -0.9;
            param.lambda = 0.0;

            // Option settings
            OpSet settings = new OpSet();
            settings.S = 101.52;
            settings.K = 100.0;
            settings.T = 0.15;
            settings.r = 0.02;
            settings.q = 0.0;
            settings.PutCall = "C";
            settings.trap = 1;

            // The Heston price
            HestonPrice HP = new HestonPrice();
            double Price = HP.HestonPriceGaussLaguerre(param,settings,x,w);
            Console.WriteLine("Heston price using 32-point Gauss Laguerre");
            Console.WriteLine("------------------------------------------ ");
            Console.WriteLine("Option Flavor =  {0,0:F5}",settings.PutCall);
            Console.WriteLine("Strike Price  =  {0,0:0}" ,settings.K);
            Console.WriteLine("Maturity      =  {0,0:F2}",settings.T);
            Console.WriteLine("Price         =  {0,0:F4}",Price);
            Console.WriteLine("------------------------------------------ ");
            Console.WriteLine(" ");
        }
    }
}
