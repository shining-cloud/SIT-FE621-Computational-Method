using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Bivariate_Characteristic_Function
{
    class BivariateCF
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
            HParam param = new HParam();
            param.kappa = 0.2;
            param.theta = 0.05;
            param.sigma = 0.3;
            param.v0 = 0.22;
            param.rho = -0.8;
            param.lambda = 0.0;

            OpSet settings = new OpSet();
            settings.S = 100.0;
            settings.T = 0.25;
            settings.r = 0.03;
            settings.q = 0.07;
            settings.PutCall = "C";
            settings.trap = 1;

            // Calculate the option price, showing identical prices under univariate and bivariate characteristic functions
            double[] K = new double[] { 90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0,101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0 };
            int N = K.Length;
            double[] SinglePriceOriginal = new double[N];
            double[] DoublePriceOriginal = new double[N];
            double[] SinglePriceTrap = new double[N];
            double[] DoublePriceTrap = new double[N];

            HestonPrice HP = new HestonPrice();
            Console.WriteLine("------------------------------------------------------");
            Console.WriteLine("          Original Heston CF      Little Trap CF");
            Console.WriteLine(" Strike  Univariate Bivariate  Univariate Bivariate");
            Console.WriteLine("------------------------------------------------------");
            for(int j=0;j<=N-1;j++)
            {
                settings.K = K[j];
                settings.trap = 0;
                SinglePriceOriginal[j] = HP.HestonPriceGaussLaguerre(param,settings,X,W,1);
                DoublePriceOriginal[j] = HP.HestonPriceGaussLaguerre(param,settings,X,W,2);
                settings.trap = 1;
                SinglePriceTrap[j] = HP.HestonPriceGaussLaguerre(param,settings,X,W,1);
                DoublePriceTrap[j] = HP.HestonPriceGaussLaguerre(param,settings,X,W,2);
                Console.WriteLine("{0,5:0} {1,10:F4} {2,10:F4} {3,10:F4} {4,10:F4}",K[j],
                    SinglePriceOriginal[j],DoublePriceOriginal[j],
                    SinglePriceTrap[j],DoublePriceTrap[j]);
            }
            Console.WriteLine("------------------------------------------------------");
        }
    }
}
