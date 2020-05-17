using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Quadratures
{
    class NewtonCotesIntegration
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] xGLa = new Double[32];
            double[] wGLa = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    xGLa[k] = double.Parse(bits[0]);
                    wGLa[k] = double.Parse(bits[1]);
                }
            }
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
            // 32-point Gauss-Lobatto Abscissas and weights
            double[] xGLo = new Double[32];
            double[] wGLo = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLobatto32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    xGLo[k] = double.Parse(bits[0]);
                    wGLo[k] = double.Parse(bits[1]);
                }
            }

            // Heston parameters
            HParam param = new HParam();
            param.kappa = 2.0;
            param.theta = 0.05;
            param.sigma = 0.3;
            param.v0 = 0.035;
            param.rho = -0.9;
            param.lambda = 0.0;

            // Option price settings
            OpSet settings = new OpSet();
            settings.S = 100.0;
            settings.T = 0.25;
            settings.r = 0.05;
            settings.q = 0.02;
            settings.trap = 1;

            // Lower and upper integration limits
            double a = 0.0;         // Lower Limit for Gauss Legendre and Newton Coates
            double b = 100.0;       // Upper Limit for Gauss Legendre
            double A = 1e-5;        // Lower Limit for Gauss Lobatto
            double B = 100.0;        // Upper Limit for Gauss Lobatto

            int N = 5000;           // Number of points for Newton Coates quadratures

            // Range of strikes and option flavor
            double[] Strikes = new double[] { 90.0,95.0,100.0,105.0,110.0 };
            settings.PutCall = "C";

            // Initialize the price vectors
            int M = Strikes.Length;
            double[] PriceGLa     = new double[M];
            double[] PriceGLe     = new double[M];
            double[] PriceGLo     = new double[M];
            double[] PriceMP      = new double[M];
            double[] PriceTrapz   = new double[M];
            double[] PriceSimp    = new double[M];
            double[] PriceSimp38  = new double[M];

            // Obtain the prices and output to console
            HestonPrice HP = new HestonPrice();
            NewtonCotes NC = new NewtonCotes();
            Console.WriteLine("Using {0:0} points for Newton Cotes and 32 points for quadrature ",N);
            Console.WriteLine(" ");
            for(int k=0;k<=M-1;k++)
            {
                settings.K = Strikes[k];
                PriceGLa[k] = HP.HestonPriceGaussLaguerre(param,settings,xGLa,wGLa);           // Gauss Laguerre
                PriceGLe[k] = HP.HestonPriceGaussLegendre(param,settings,xGLe,wGLe,a,b);       // Gauss Legendre
                PriceGLo[k] = HP.HestonPriceGaussLegendre(param,settings,xGLo,wGLo,A,B);       // Gauss Lobatto
                PriceMP[k]     = NC.HestonPriceNewtonCotes(param,settings,1,a,b,N);           // Mid point rule
                PriceTrapz[k]  = NC.HestonPriceNewtonCotes(param,settings,2,A,B,N);           // Trapezoidal rule
                PriceSimp[k]   = NC.HestonPriceNewtonCotes(param,settings,3,A,B,N);           // Simpson's Rule
                PriceSimp38[k] = NC.HestonPriceNewtonCotes(param,settings,4,A,B,N);           // Simpson's 3/8 rule
                Console.Write("Strike {0,2}", Strikes[k]);
                Console.WriteLine(" ------------------");
                Console.WriteLine("Mid Point       {0:F5}",PriceMP[k]);
                Console.WriteLine("Trapezoidal     {0:F5}",PriceTrapz[k]);
                Console.WriteLine("Simpson's       {0:F5}",PriceSimp[k]);
                Console.WriteLine("Simpson's 3/8   {0:F5}",PriceSimp38[k]);
                Console.WriteLine("Gauss Laguerre  {0:F5}",PriceGLa[k]);
                Console.WriteLine("Gauss Legendre  {0:F5}",PriceGLe[k]);
                Console.WriteLine("Gauss Lobatto   {0:F5}",PriceGLo[k]);
                Console.WriteLine(" ");
            }
        }
    }
}


/*
namespace Heston_Quadratures
{
    public partial class Quadratures
    {
        static void Main(string[] args)
        {

            // Heston parameters
            HParam param = new HParam();
            param.kappa = 5.0;
            param.theta = 0.16;
            param.sigma = 0.9;
            param.v0 = 0.0625;
            param.rho = 0.1;
            param.lambda = 0.0;

            // Clarke and Parrott option price settings, but for European puts
            OpSet settings = new OpSet();
            settings.K = 10.0;
            settings.T = 0.25;
            settings.r = 0.1;
            settings.q = 0.0;
            settings.trap = 1;

            // Lower and upper integration limits
            double A = 1e-10;        // Lower Limit
            double B = 200.0;        // Upper Limit 
            int N = 50000;           // Number of points 

            // Range of strikes and option flavor
            double[] Spot= new double[] { 8.0,9.0,10.0,11.0,12.0 };
            settings.PutCall = "P";

            // Initialize the price vectors
            int M = Spot.Length;
            double[] PriceTrapz  = new double[M];
            double[] PriceSimp   = new double[M];
            double[] PriceSimp38 = new double[M];

            // Obtain the prices and output to console
            for(int k=0;k<=M-1;k++)
            {
                settings.S = Spot[k];
                PriceTrapz[k]  = HestonPriceNewtonCoates(param,settings,2,A,B,N);           // Trapezoidal rule
                PriceSimp[k]   = HestonPriceNewtonCoates(param,settings,3,A,B,N);           // Simpson's Rule
                PriceSimp38[k] = HestonPriceNewtonCoates(param,settings,4,A,B,N);           // Simpson's 3/8 rule
                Console.Write("Strike {0,2}",Spot[k]);
                Console.WriteLine(" ------------------");
                Console.WriteLine("Trapezoidal     {0:F6}",PriceTrapz[k]);
                Console.WriteLine("Simpson's       {0:F6}",PriceSimp[k]);
                Console.WriteLine("Simpson's 3/8   {0:F6}",PriceSimp38[k]);
                Console.WriteLine(" ");
            }
        }
    }
}

*/

