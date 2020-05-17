using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Fundamental_Transform
{
    class LewisFT
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
            double K = 100.0;				    // Strike Price
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

            // Heston Price
            HestonPrice HP = new HestonPrice();
            double HestonPrice = HP.HestonPriceGaussLaguerre(PutCall,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);

            // Lewis Price form 1
            LewisPrice LP = new LewisPrice();
            double ki1 = 1.5;
            int form1 = 1;
            double LewisPrice1 = LP.HestonLewisCallPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w);

            // Lewis Price form 2
            double ki2 = 0.5;
            int form2 = 2;
            double LewisPrice2 = LP.HestonLewisCallPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w);

            // Write the results
            Console.WriteLine("Lewis Fundamental Transform");
            Console.WriteLine("Method           Price   ");
            Console.WriteLine("---------------------------");
            Console.WriteLine("Heston           {0:F4}",HestonPrice);
            Console.WriteLine("Lewis form 1     {0:F4}",LewisPrice1);
            Console.WriteLine("Lewis form 2     {0:F4}",LewisPrice2);
            Console.WriteLine("---------------------------");
        }
    }
}


