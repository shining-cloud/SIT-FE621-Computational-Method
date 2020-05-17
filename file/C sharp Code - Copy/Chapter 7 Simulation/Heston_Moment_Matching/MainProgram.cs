using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Moment_Matching
{
    class HestonMomentMatching
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
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
            HParam param = new HParam();
            param.kappa = 6.2;
            param.theta = 0.06;
            param.sigma = 0.5;
            param.v0 = 0.03;
            param.rho = -0.7;
            param.lambda = 0.0;

            OpSet settings = new OpSet();
            settings.S = 100.0;
            settings.K = 90.0;
            settings.T = 0.25;
            settings.r = 0.03;
            settings.q = 0.02;
            settings.PutCall = "C";
            settings.trap = 1;

            // Simulation settings
            int NT = 100;
            int NS = 5000;
            
            // Calculate the IJK price
            MMSimulation MM = new MMSimulation();
            double MomentPrice = MM.MMPrice(param,settings,NT,NS);

            // Calculate the closed-form European option price
            HestonPrice HP = new HestonPrice();
            double ClosedPrice = HP.HestonPriceGaussLaguerre(param,settings,x,w);

            // Calculate the errors;
            double MMError = Math.Abs((MomentPrice-ClosedPrice)/ClosedPrice*100);

            // Output the results
            Console.WriteLine("European Heston price by simulation --------------------");
            Console.WriteLine("Uses  {0:0} stock price paths, and {1:0} time increments",NS,NT);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Closed form price        {0:F5}",ClosedPrice);
            Console.WriteLine("Moment Matching price    {0:F5}   Percent Error {1:F5}",MomentPrice,MMError);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine(" ");
        }
    }
}
