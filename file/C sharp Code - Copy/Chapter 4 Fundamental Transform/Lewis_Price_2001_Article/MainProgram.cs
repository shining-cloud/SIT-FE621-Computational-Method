using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Price_2001_Article
{
    class Lewis2001
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
            double S = 100.0;				    // Spot Price
            double T = 0.25;			        // Maturity in Years
            double r = 0.05;					// Interest Rate
            double q = 0.01;                    // Dividend yield
            double rho = -0.9;					// Heston Parameter: Correlation
            double kappa = 2;					// Heston Parameter 
            double theta = 0.05;				// Heston Parameter 
            double lambda = 0.0;				// Heston Parameter 
            double sigma = 0.1;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.05;					// Heston Parameter: Current Variance
            string PutCall = "C";               // "P"ut or "Call"
            int trap = 1;                       // 1="Little Trap" characteristic function

            // Calculate the Heston and Lewis Prices for a range of strikes
            HestonPrice HP = new HestonPrice();
            LewisPrice LP = new LewisPrice();
            double[] K = new double[11] { 95.0,96.0,97.0,98.0,99.0,100.0,101.0,102.0,103.0,104.0,105.0 };
            double[] HestonPrice = new Double[11];
            double[] LewisPrice  = new Double[11];
            double[] Error       = new Double[11];
            Console.WriteLine("Lewis Price of Heston model by Fundamental transform");
            Console.WriteLine("-------------------------------------------------------");
            Console.WriteLine("Strike   Heston Price    Lewis Price    Percent Error");
            Console.WriteLine("-------------------------------------------------------");
            for(int j=0;j<=10;j++)
            {
                HestonPrice[j] = HP.HestonPriceGaussLaguerre(PutCall,S,K[j],r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
                LewisPrice[j] = LP.LewisPrice311(S,K[j],r,q,T,theta,kappa,sigma,rho,v0,trap,x,w);
                Error[j] = (LewisPrice[j] - HestonPrice[j]) / HestonPrice[j] * 100.0;
                Console.WriteLine("{0,4} {1,12:F4} {2,15:F4} {3,15:F2} ", K[j], HestonPrice[j], LewisPrice[j], Error[j]);
            }
            Console.WriteLine("-------------------------------------------------------");
        }
    }
}
