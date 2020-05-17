using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Transformed_Volatility
{
    public partial class TV
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
            settings.T = 3.0/12.0;
            settings.r = 0.03;
            settings.q = 0.02;
            settings.PutCall = "C";
            settings.trap = 1;

            // Simulation settings
            int NT = 100;
            int NS = 2000;

            // Calculate the Euler price and the transformed volatility price
            TVSimulation TV = new TVSimulation();
            double TVolPrice  = TV.TransVolPrice("TV",param,settings,NT,NS);
            double EulerPrice = TV.TransVolPrice("Euler",param,settings,NT,NS);

            // Calculate the closed-form European option price
            double ClosedPrice = HestonPriceGaussLaguerre(param,settings,x,w);

            // Calculate the errors;
            double TVolError  = Math.Abs((TVolPrice-ClosedPrice)/ClosedPrice*100);
            double EulerError = Math.Abs((EulerPrice-ClosedPrice)/ClosedPrice*100);

            // Output the results
            Console.WriteLine("Uses  {0:0} stock price paths and {1:0} time increments",NS,NT);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Method                       Price     PercentError");
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Exact with Gauss Laguerre  {0,10:F5}",ClosedPrice);
            Console.WriteLine("Euler discretization       {0,10:F5}  {1,10:F5}",EulerPrice,EulerError);
            Console.WriteLine("Transformed volatility     {0,10:F5}  {1,10:F5}",TVolPrice,TVolError);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine(" ");
        }
    }
}
