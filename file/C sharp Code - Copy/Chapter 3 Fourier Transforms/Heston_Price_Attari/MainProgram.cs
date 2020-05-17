using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Price_Attari
{
    class HestonAttariPrice
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
            double S = 30.0;				    // Spot Price
            double K = 20.0;				    // Strike Price
            double T = 1.0/12.0;		        // Maturity in Years
            double r = 0.01;					// Interest Rate
            double q = 0.0;                     // Dividend yield
            double kappa = 1.4;					// Heston Parameter 
            double theta = 0.05;				// Heston Parameter 
            double sigma = 0.3;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.05;					// Heston Parameter: Current Variance
            double rho = -0.8;					// Heston Parameter: Correlation
            double lambda = 0.0;				// Heston Parameter 
            string PutCall = "C";               // "P"ut or "Call"
            int trap = 1;                       // 1="Little Trap" characteristic function

            // Calculate the Heston and Attari prices
            HestonAttari HA = new HestonAttari();
            double HestonPrice = HA.HestonPriceGaussLaguerre(PutCall,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double AttariPrice = HA.AttariPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);

            // Output the results
            Console.WriteLine("Method      Call Price");
            Console.WriteLine("----------------------");
            Console.WriteLine("Heston      {0:F8}",HestonPrice);
            Console.WriteLine("Attari      {0:F8}",AttariPrice);
            Console.WriteLine("----------------------");
        }
    }
}


