using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Euler_and_Milstein
{
    class EulerMilstein
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
            param.kappa = 2.0;
            param.theta = 0.05;
            param.sigma = 0.3;
            param.v0 = 0.035;
            param.rho = -0.9;
            param.lambda = 0.0;

            OpSet settings = new OpSet();
            settings.S = 100.0;
            settings.K = 100.0;
            settings.T = 0.1;
            settings.r = 0.14;
            settings.q = 0.10;
            settings.PutCall = "P";
            settings.trap = 1;

            // Simulation settings
            string negvar = "Truncation";
            double alpha = 0.5;
            int NT = 100;
            int NS = 5000;
            
            // Calculate the Euler price
            Simulation SI = new Simulation();
            HestonPrice HP = new HestonPrice();
            double EulerPrice = SI.EulerMilsteinPrice("Euler",negvar,param,settings,alpha,NT,NS,settings.PutCall);

            // Calculate the Milstein price
            double MilsteinPrice = SI.EulerMilsteinPrice("Milstein",negvar,param,settings,alpha,NT,NS,settings.PutCall);

            // Calculate the Implicit Milstein price
            double IMPrice = SI.EulerMilsteinPrice("IM",negvar,param,settings,alpha,NT,NS,settings.PutCall);

            // Calculate the Weighted Implicit Explicit Milstein price
            double WMPrice = SI.EulerMilsteinPrice("WM",negvar,param,settings,alpha,NT,NS,settings.PutCall);

            // Calculate the closed-form European option price
            double ClosedPrice = HP.HestonPriceGaussLaguerre(param,settings,x,w);

            // Calculate the errors;
            double EulerError    = Math.Abs((EulerPrice   -ClosedPrice)/ClosedPrice*100);
            double MilsteinError = Math.Abs((MilsteinPrice-ClosedPrice)/ClosedPrice*100);
            double IMError       = Math.Abs((IMPrice      -ClosedPrice)/ClosedPrice*100);
            double WMError       = Math.Abs((WMPrice      -ClosedPrice)/ClosedPrice*100);

            // Output the results
            Console.WriteLine("European Heston price by simulation --------------------");
            Console.WriteLine("Uses {0:0} stock price paths and {1:0} time increments", NS,NT);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("                         Value   PercentError");
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Closed form price        {0:F5}",ClosedPrice);
            Console.WriteLine("Euler price              {0:F5}  {1,8:F5}",EulerPrice,EulerError);
            Console.WriteLine("Milstein price           {0:F5}  {1,8:F5}",MilsteinPrice,MilsteinError);
            Console.WriteLine("Implicit Milstein price  {0:F5}  {1,8:F5}",IMPrice,IMError);
            Console.WriteLine("Weighted Milstein price  {0:F5}  {1,8:F5}",WMPrice,WMError);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine(" ");
        }
        // Random number in (a,b) ==========================================================
        private static readonly Random U = new Random();
        private static readonly object sync = new object();
        public static double RandomNum(double a,double b)
        {
            int divisor = 1000000000;
            lock(sync) { return a + (b-a)*U.Next(0,divisor)/divisor; }
        }

    }
}
