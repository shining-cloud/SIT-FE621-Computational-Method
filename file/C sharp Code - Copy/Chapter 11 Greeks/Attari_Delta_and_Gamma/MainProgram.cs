using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Attari_Delta_and_Gamma
{
    class AttariDeltaGamma
    {
        static void Main(string[] args)
        {
            AttariGreeks AG = new AttariGreeks();
            AttariPrice AP = new AttariPrice();

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
            double S = 25.0;				    // Spot Price
            double T = 0.5;		                // Maturity in Years
            double r = 0.05;					// Interest Rate
            double q = 0.0;                     // Dividend yield
            double kappa = 1.4;					// Heston Parameter 
            double theta = 0.05;				// Heston Parameter 
            double sigma = 0.5;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.01;					// Heston Parameter: Current Variance
            double rho = -0.8;					// Heston Parameter: Correlation
            double lambda = 0.0;				// Heston Parameter 
            string PutCall = "C";               // "P"ut or "Call"
            int trap = 1;                       // 1="Little Trap" characteristic function

            // Define the range of strikes
            double[] K = new double[] { 20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0 };
            int NK = K.Length;
            double[] DeltaFD = new double[NK];
            double[] GammaFD = new double[NK];
            double[] Delta = new double[NK];
            double[] Gamma = new double[NK];
            double dS = 0.1;

            // Calculate delta and gamma using closed form and finite differences, and write to screen
            double Price0,Price1,Price1_;
            Console.WriteLine("---------------------------------------------------");
            Console.WriteLine("Strike   DeltaFD     Delta      GammaFD     Gamma");
            Console.WriteLine("---------------------------------------------------");
            for(int k=0;k<=NK-1;k++)
            {
                Price1  = AP.AttariPriceGaussLaguerre(PutCall,S+dS,K[k],T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
                Price0  = AP.AttariPriceGaussLaguerre(PutCall,S   ,K[k],T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
                Price1_ = AP.AttariPriceGaussLaguerre(PutCall,S-dS,K[k],T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
                DeltaFD[k] = (Price1 - Price1_)/(2.0*dS);
                GammaFD[k] = (Price1 - 2*Price0 + Price1_)/(dS*dS);
                Delta[k] = AG.AttariGreeksAnalytic(PutCall,kappa,theta,lambda,rho,sigma,T,K[k],S,r,q,v0,trap,"Delta",x,w);
                Gamma[k] = AG.AttariGreeksAnalytic(PutCall,kappa,theta,lambda,rho,sigma,T,K[k],S,r,q,v0,trap,"Gamma",x,w);
                Console.WriteLine("{0,5:F0} {1,10:F4} {2,10:F4} {3,10:F4} {4,10:F4}",K[k],DeltaFD[k],Delta[k],GammaFD[k],Gamma[k]);
            }
            Console.WriteLine("---------------------------------------------------");
        }
    }
}
