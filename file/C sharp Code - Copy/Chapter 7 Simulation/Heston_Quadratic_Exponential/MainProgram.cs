using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Quadratic_Exponential
{
    public partial class QE
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
            settings.K = 100.0;
            settings.T = 0.25;
            settings.r = 0.03;
            settings.q = 0.02;
            settings.PutCall = "C";
            settings.trap = 1;

            // Simulation settings
            double gamma1 = 0.5;
            double gamma2 = 0.5;
            double phic = 1.5;
            int MC = 1;
            int NT = 100;
            int NS = 1500;
            
            // Calculate the IJK price
            QESimulation QE = new QESimulation();
            double QuadExpPrice = QE.QEPrice(param,settings,gamma1,gamma2,NT,NS,MC,phic,settings.PutCall);

            // Calculate the closed-form European option price
            HestonPrice HP = new HestonPrice();
            double ClosedPrice = HP.HestonPriceGaussLaguerre(param,settings,x,w);

            // Calculate the errors;
            double QError = Math.Abs((QuadExpPrice-ClosedPrice)/ClosedPrice*100);

            // Output the results
            Console.WriteLine("European Heston price by simulation --------------------");
            Console.WriteLine("Uses  {0:0} stock price paths and {1:0} time increments",NS,NT);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Closed form price        {0:F5}",ClosedPrice);
            Console.WriteLine("Quadratic exp price      {0:F5}      Error {1,5:F5}",QuadExpPrice,QError);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine(" ");
        }
    }
}
