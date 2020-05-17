using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace MultiDomain_Integration_Zhu
{
    class MultiDomainZhu
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Legendre Abscissas and weights
            double[] xGLe = new Double[32];
            double[] wGLe = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLegendre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    xGLe[k] = double.Parse(bits[0]);
                    wGLe[k] = double.Parse(bits[1]);
                }
            }

            // Heston parameters
            HParam param = new HParam();
            param.kappa = 2.0;
            param.theta = 0.05;
            param.sigma = 0.3;
            param.v0 = 0.05;
            param.rho = 0.45;
            param.lambda = 0.0;

            // Option price settings
            OpSet settings = new OpSet();
            settings.S = 100.0;
            settings.K = 100.0;
            settings.r = 0.05;
            settings.q = 0.01;
            settings.trap = 1;
            settings.PutCall = "C";
            settings.T = 0.5;

            // Price by Simpson's rule
            NewtonCotesPrice NCP = new NewtonCotesPrice();
            int method = 3;
            double a = 1e-10;
            double b = 150;
            int NS = 10000;
            double PriceSimpson = NCP.HestonPriceNewtonCotes(param,settings,method,a,b,NS);

            // Potential integration domain 
            double lo = 1e-10;
            double hi = 150;
            int N = 1000;
            double dA = (hi-lo)/Convert.ToDouble(N);
            double[] A = new double[N];
            for(int j=0;j<=N-1;j++)
                A[j] = lo + j*dA;

            // The modified domain Gauss-Legendre price
            HestonPriceMD HPMD = new HestonPriceMD();
            double tol = 1e-6;
            OutputMD output = HPMD.HestonPriceGaussLegendreMD(param,settings,xGLe,wGLe,A,tol);
            double PriceMD = output.Price;
            double lower = output.lower;
            double upper = output.upper;
            int Npoints = output.Npoints;
            double errorMD = PriceMD - PriceSimpson;

            // The price using Newton Cotes
            double PriceNC = NCP.HestonPriceNewtonCotes(param,settings,method,lower,upper,Npoints);
            double errorNC = PriceNC - PriceSimpson;

            // Write the results
            Console.WriteLine("Zhu Multi-Domain Integration Algorithm");
            Console.WriteLine("-------------------------------------------");
            Console.WriteLine("Multi-Domain tolerance      {0:E3}",tol);
            Console.WriteLine("Lower integration limit     {0:E3}",lower);
            Console.WriteLine("Upper integration limit     {0:F3}",upper);
            Console.WriteLine("Number of integration pts   {0:0}",Npoints);
            Console.WriteLine("-------------------------------------------");
            Console.WriteLine("Method                 Price        Error ");
            Console.WriteLine("-------------------------------------------");
            Console.WriteLine("{0:0}-point Simpson  {1,10:F4}",N,PriceSimpson);
            Console.WriteLine("Multi-Domain        {0,10:F4} {1,10:F4}",PriceMD,errorMD);
            Console.WriteLine("Newton-Cotes        {0,10:F4} {1,10:F4}",PriceNC,errorNC);
            Console.WriteLine("-------------------------------------------");
        }
    }
}
