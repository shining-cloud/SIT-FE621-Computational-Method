using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Double_Heston_Gauthier_Possamai
{
    class DoubleHestonGauthierPossamai
    {
        static void Main(string[] args)
        {
            HestonPriceDH HPDH = new HestonPriceDH();

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
            // Double Heston parameters
            DHParam param = new DHParam();
            param.v01 = 0.6*0.6;
            param.v02 = 0.7*0.7;
            param.sigma1 = 0.10;
            param.sigma2 = 0.20;
            param.kappa1 = 0.90;
            param.kappa2 = 1.20;
            param.rho1 = -0.5;
            param.rho2 = -0.5;
            param.theta1 = 0.10;
            param.theta2 = 0.15;

            // Option settings
            OpSet settings = new OpSet();
            settings.S = 61.90;
            settings.K = 100.0;
            settings.r = 0.03;
            settings.q = 0.0;
            settings.PutCall = "C";
            settings.trap = 1;

            // Strikes and Maturities
            double[] percent = new double[6] { 1.0,1.0,0.7,0.7,1.3,1.3 };
            double[] Mat     = new double[6] { 1.0,10.0,1.0,10.0,1.0,10.0 };

            // Calculate the option prices using Gauss Laguerre integration
            double[] Christoffersen = new double[6];
            double[] Gauthier       = new double[6];

            Console.WriteLine("-------------------------------------------");
            Console.WriteLine("Table 3 from Gauthier and Possamai");
            Console.WriteLine(" ");
            Console.WriteLine("Strike  Maturity   Christoffersen  Gauthier");
            Console.WriteLine("-------------------------------------------");
            for(int k=0;k<=5;k++)
            {
                settings.T = Mat[k];
                settings.K = settings.S * percent[k];
                settings.trap = 0;
                Christoffersen[k] = HPDH.DoubleHestonPriceGaussLaguerre(param,settings,X,W);
                settings.trap = 1;
                Gauthier[k]       = HPDH.DoubleHestonPriceGaussLaguerre(param,settings,X,W);
                Console.WriteLine("{0,3:0.00} {1,5:0} {2,15:F5} {3,15:F5}",settings.K, settings.T,Christoffersen[k],Gauthier[k]);
            }
            Console.WriteLine("-------------------------------------------");
        }
    }
}
