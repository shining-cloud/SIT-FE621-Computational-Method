using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_FFT
{
    class MainProgram
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
            // Heston parameters
            HParam param = new HParam();
            param.kappa = 0.2;
            param.theta = 0.05;
            param.sigma = 0.3;
            param.v0 = 0.05;
            param.rho = -0.7;
            param.lambda = 0.0;

            // Option settings
            OpSet settings = new OpSet();
            settings.S = 50.0;
            settings.T = 0.5;
            settings.r = 0.03;
            settings.q = 0.05;
            settings.PutCall = "C";
            settings.trap = 1;

            // FFT settings
            double M = Math.Pow(2.0,10.0);
            int N = Convert.ToInt16(M);
            double alpha = 1.5;
            double uplimit = 100.0;
            string rule = "Trapezoidal";
            double[,] PriceFFT = new double[N,2];

            // Discount the spot price by the dividend yield, use as input in FFT
            settings.S = settings.S*Math.Exp(-settings.q*settings.T);

            // Obtain the FFT strikes and calls
            FFT FFT = new FFT();
            PriceFFT = FFT.HestonFFT(N,uplimit,alpha,rule,param,settings);

            // Obtain the extact price and compare to FFT price
            HestonPrice HP = new HestonPrice();
            double[] PriceHeston = new double[N];
            for(int j=0;j<=N-1;j++)
            {
                settings.K = PriceFFT[j,0];
                PriceHeston[j] = HP.HestonPriceGaussLaguerre(param,settings,x,w);
            }

            // Output the results near the money
            int start = N/2 - 5;
            int end   = N/2 + 5;
            Console.WriteLine("Fast Fourier Transform");
            Console.WriteLine("----------------------------------------");
            Console.WriteLine("Strike        FFT Price     Heston Price");
            Console.WriteLine("----------------------------------------");
            for(int j=start;j<=end;j++)
                Console.WriteLine("{0,2:F2} {1,15:F4} {2,15:F4} ",PriceFFT[j,0],PriceFFT[j,1],PriceHeston[j]);

            Console.WriteLine("----------------------------------------");
        }
    }
}


