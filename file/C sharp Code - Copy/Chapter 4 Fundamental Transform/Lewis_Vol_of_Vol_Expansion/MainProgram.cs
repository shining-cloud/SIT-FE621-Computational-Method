using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Vol_of_Variance_Expansion
{
    class LewisVolVar
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
            double rf = 0.0;					// Interest Rate
            double q = 0.0;                     // Dividend yield

            double rho = -0.5;					// Heston Parameter: Correlation
            double kappa = 4.0; 			    // Heston Parameter 
            double theta = 0.09/4.0;			// Heston Parameter 
            double lambda = 0.0;				// Heston Parameter 
            double sigma = 0.1;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.0225;					// Heston Parameter: Current Variance
            string PutCall = "C";               // "P"ut or "Call"
            int trap = 1;                       // 1="Little Trap" characteristic function
            int NK = 7;
            double[] K = new double[7]          // Strikes
            { 70.0,80.0,90.0,100.0,110.0,120.0,130.0 };

            // Time average of the deterministic variance and volatility
            double v = theta + (v0-theta)*(1.0-Math.Exp(-kappa*T))/(kappa*T);
            double IV = Math.Sqrt(v);

            // Classes
            BlackScholes BS = new BlackScholes();
            HestonPrice HP = new HestonPrice();
            Lewis L = new Lewis();
            BisectionAlgo BA = new BisectionAlgo();

            // Table 3.3.1 on Page 81 of Lewis (2001)
            double[] SeriesIPrice = new double[NK];         // Series I price and implied vol
            double[] IV1 = new double[NK];
            double[] SeriesIIPrice = new double[NK];        // Series II price and implied vol
            double[] IV2 = new double[NK];
            double[] ExactPrice = new double[NK];           // Exact Heston price and implied vol
            double[] IVe = new double[NK]; 
            double[] BSPrice = new double[NK];              // Black Scholes price

            // Bisection algorithm settings
            double a = 0.001;
            double b = 5.0;
            double Tol = 1.0e-5;
            int MaxIter = 1000;

            for(int k=0;k<=NK-1;k++)
            {
                SeriesIPrice[k] = L.SeriesICall(S,K[k],rf,q,T,v0,rho,theta,kappa,sigma);
                IV1[k] = BA.BisecBSIV(PutCall,S,K[k],rf,q,T,a,b,SeriesIPrice[k],Tol,MaxIter);
                double[] SeriesII = L.SeriesIICall(S,K[k],rf,q,T,v0,rho,theta,kappa,sigma);
                SeriesIIPrice[k] = SeriesII[0];
                IV2[k] = SeriesII[1];
                ExactPrice[k] = HP.HestonPriceGaussLaguerre(PutCall,S,K[k],rf,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
                BSPrice[k] = BS.BSC(S,K[k],rf,q,v,T);
            }
            Console.WriteLine("Lewis Vol of Vol expansion");
            Console.WriteLine("-----------------------------------------------");
            Console.WriteLine("Strike Price         70       80       90    100    110    120    130");
            Console.WriteLine("Exact {0,20:F4} {1,8:F4} {2,8:F4}",ExactPrice[0],ExactPrice[1],ExactPrice[2]);




        }
    }
}
