using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Alfonsi
{
    class Program
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] X = new Double[32];
            double[] W = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    X[k] = double.Parse(bits[0]);
                    W[k] = double.Parse(bits[1]);
                }
            }

            // Option features
            double r = 0.03;           // Risk free rate
            double q = 0.02;           // Dividend yield
            double Mat = .25;          // Maturity in years
            double S0 = 100;           // Spot price
            double K  = 90;            // Strike price
            string PutCall = "C";      // 'P'ut or 'C'all
            OpSet opset = new OpSet();
            opset.S = S0;
            opset.K = K;
            opset.r = r;
            opset.q = q;
            opset.T = Mat;
            opset.PutCall = PutCall;
            opset.trap = 1;

            // Heston parameters
            HParam param = new HParam();
            param.kappa =  1.2;      // Variance reversion speed
            param.theta =  0.06;     // Variance reversion level
            param.sigma =  0.5;      // Volatility of Variance
            param.rho   = -0.70;     // Correlation between Brownian motions
            param.v0    =  0.03;     // Initial variance
            param.lambda = 0.0;      // Risk

            // Simulation features
            int N = 5000;             // Number of stock price paths
            int T = 100;              // Number of time steps per path

            // Closed form Heston
            HestonAnalytics HPrice = new HestonAnalytics();
            double HesPrice = HPrice.HestonPriceGaussLaguerre(param,opset,X,W);
            
            // Alfonsi Simulated price
            AlfonsiSim APrice = new AlfonsiSim();
            double SimPrice = APrice.AlfonsiPrice(param,S0,K,Mat,r,q,T,N,PutCall);

            // Error
            double DollarError = HesPrice - SimPrice;

            // Output the results
            Console.WriteLine("Alfonsi simulated price using {0:0} simulations and {1:0} time steps",N,T);
            Console.WriteLine("-----------------------------------------------");
            Console.WriteLine("Method          Price            Dollar Error");
            Console.WriteLine("-----------------------------------------------");
            Console.WriteLine("Closed Form     {0:F5}",HesPrice);
            Console.WriteLine("Alfonsi         {0:F5}        {1:F5}",SimPrice,DollarError);
            Console.WriteLine("-----------------------------------------------");
        }
    }
}
