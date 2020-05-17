using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Kahl_Jackel_Integration
{
    class KahlJackelIntegration
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
            param.kappa = 1.0;
            param.theta = 0.06;
            param.sigma = 0.5;
            param.v0 = 0.06;
            param.rho = -0.8;
            param.lambda = 0.0;

            // Option price settings
            OpSet settings = new OpSet();
            settings.S = 10.0;
            settings.K = 7.0;
            settings.T = 1.0/12.0;
            settings.r = 0.06;
            settings.q = 0.04;
            settings.trap = 1;
            settings.PutCall = "C";

            // Lower and upper integration limits
            double a = 0.0;         // Lower Limit for Gauss Legendre
            double b = 100.0;       // Upper Limit for Gauss Legendre
            double A = 1e-5;        // Lower Limit for Gauss Lobatto
            double B = 100.0;       // Upper Limit for Gauss Lobatto

            // Initialize the price vectors
            HestonPrice HP = new HestonPrice();
            KahlJackel KJ = new KahlJackel();
            double PriceGLa = HP.HestonPriceGaussLaguerre(param,settings,xGLa,wGLa);           // Heston Gauss Laguerre
            double PriceGLe = HP.HestonPriceGaussLegendre(param,settings,xGLe,wGLe,a,b);       // Heston Gauss Legendre
            double PriceGLo = HP.HestonPriceGaussLegendre(param,settings,xGLo,wGLo,A,B);       // Heston Gauss Lobatto
            double PriceKJ  = KJ.HestonPriceKahlJackel(param,settings,xGLo,wGLo);              // Kahl Jackel Gauss Lobatto

            Console.WriteLine("Kahl Jackel Integration scheme");
            Console.WriteLine("Method                       Price");
            Console.WriteLine("-------------------------------------");
            Console.WriteLine("Heston Gauss Laguerre       {0:F5}",PriceGLa);
            Console.WriteLine("Heston Gauss Legendre       {0:F5}",PriceGLe);
            Console.WriteLine("Heston Gauss Lobatto        {0:F5}",PriceGLo);
            Console.WriteLine("Kahl Jackel Gauss Lobatto   {0:F5}",PriceKJ);
            Console.WriteLine("-------------------------------------");
        }
    }
}

