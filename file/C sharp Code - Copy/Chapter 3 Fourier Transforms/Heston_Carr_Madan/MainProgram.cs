using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Carr_Madan
{
    class HestonCM
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
            double T = 0.5 ;			        // Maturity in Years
            double r = 0.10;					// Interest Rate
            double q = 0.07;                    // Dividend yield
            double rho = -0.7;					// Heston Parameter: Correlation
            double kappa = 2.0;					// Heston Parameter 
            double theta = 0.06;				// Heston Parameter 
            double lambda = 0.0;				// Heston Parameter 
            double sigma = 0.1;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.06;					// Heston Parameter: Current Variance
            int trap = 1;                       // 1="Little Trap" characteristic function
            double alpha = 1.75;                // Carr Madan damping factor

            // Calculate the prices using the Heston and Carr-Madan integrands
            HestonPrice HP = new HestonPrice();
            double HestonCall    = HP.HestonPriceGaussLaguerre("Heston","C",alpha,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double HestonPut     = HP.HestonPriceGaussLaguerre("Heston","P",alpha,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double CarrMadanCall = HP.HestonPriceGaussLaguerre("CarrMadan","C",alpha,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double CarrMadanPut  = HP.HestonPriceGaussLaguerre("CarrMadan","P",alpha,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);

            // Output the results
            Console.WriteLine(" ");
            Console.WriteLine("Method              Call Price");
            Console.WriteLine("-----------------------------");
            Console.WriteLine("Heston Integrand      {0:F4} ", HestonCall);
            Console.WriteLine("Carr-Madan Integrand  {0:F4} ", CarrMadanCall);
            Console.WriteLine("-----------------------------");
            Console.WriteLine(" ");
            Console.WriteLine("Method               Put Price");
            Console.WriteLine("-----------------------------");
            Console.WriteLine("Heston Integrand      {0:F4} ",HestonPut);
            Console.WriteLine("Carr-Madan Integrand  {0:F4} ",CarrMadanPut);
            Console.WriteLine("-----------------------------");
            Console.WriteLine(" ");
        }
    }
}

