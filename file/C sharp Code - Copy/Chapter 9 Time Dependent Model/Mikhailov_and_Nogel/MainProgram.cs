using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Mikhailov_and_Nogel
{
    class MikhailovNogelPrice
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
            HParam param = new HParam();
            param.kappa = 1.0;
            param.theta = 0.1;
            param.sigma = 0.2;
            param.v0 = 0.1;
            param.rho = -0.3;

            OpSet settings = new OpSet();
            settings.S = 1.0;
            settings.r = 0.0;
            settings.q = 0.0;
            settings.PutCall = "C";
            settings.trap = 1;

            // Past maturities
            double[] tau0 = new double[] { 5.0/3.0, 5.0/3.0};
            
            // Current maturity
            double tau = 5.0/3.0;

            // Past parameter values for kappa=4 and kappa=2
            double[,] param0 = { { 4.0,param.theta,param.sigma,param.v0,param.rho },
                                 { 2.0,param.theta,param.sigma,param.v0,param.rho }}; 

            // Current parameter values
            double[] K = new double[]{0.5,0.75,1.0,1.25,1.5};
            int N = K.Length;
            double[] NMPrice = new double[N];

            // Calculate the Time Dependent prices of Mikhailov and Nogel
            HestonPriceTD HPTD = new HestonPriceTD();
            for(int j=0;j<=N-1;j++)
            {
                settings.K = K[j];
                NMPrice[j] = HPTD.MNPriceGaussLaguerre(param,param0,tau,tau0,settings,x,w);
            }

            // For comparison, the time-independent (constant parameter) prices
            // Use average value of kappa and a maturity of 5 years
            param.kappa = 7.0/3.0;
            double[] tau00 = new double[] {0.0,0.0};
            double[] PriceInd = new double[N];
            tau = 5.0;
            for(int j=0;j<=N-1;j++)
            {
                settings.K = K[j];
                PriceInd[j] = HPTD.MNPriceGaussLaguerre(param,param0,tau,tau00,settings,x,w);
            }

            // Output the results
            Console.WriteLine("Strike       Mikhailov-Nogel TD Price   Static Parameters Price");
            Console.WriteLine("---------------------------------------------------------------");
            for(int j=0;j<=N-1;j++)
            {
                Console.WriteLine("{0:0.00} {1,20:F5} {2,25:F5}",K[j],NMPrice[j],PriceInd[j]);
            }
            Console.WriteLine("---------------------------------------------------------------");
        }
    }
}
