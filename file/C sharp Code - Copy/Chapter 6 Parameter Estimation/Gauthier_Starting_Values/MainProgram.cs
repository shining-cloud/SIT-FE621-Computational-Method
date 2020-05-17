using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Gauthier_Starting_Values
{
    class GauthierStartingValues
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] X = new Double[32];
            double[] W = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    X[k] = double.Parse(bits[0]);
                    W[k] = double.Parse(bits[1]);
                }

            // Size of strikes (NK) and maturities (NT)
            int NK = 7;
            int NT = 4;

            // Market prices of SPX puts
            double[,] MktPrice = new Double[NK,NT];
            using(TextReader reader = File.OpenText("../../SPX_Puts.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    MktPrice[k,0] = double.Parse(bits[0]);
                    MktPrice[k,1] = double.Parse(bits[1]);
                    MktPrice[k,2] = double.Parse(bits[2]);
                    MktPrice[k,3] = double.Parse(bits[3]);
                }

            // SPX put implied volatility
            double[,] MktIV = new Double[NK,NT];
            using(TextReader reader = File.OpenText("../../SPX_PutIV.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    MktIV[k,0] = double.Parse(bits[0]);
                    MktIV[k,1] = double.Parse(bits[1]);
                    MktIV[k,2] = double.Parse(bits[2]);
                    MktIV[k,3] = double.Parse(bits[3]);
                }

            // Option settings
            OpSet settings = new OpSet();
            settings.S = 137.14;
            settings.r = 0.00;
            settings.q = 0.00;
            settings.trap = 1;
            double[] K = new Double[7] { 120.0,125.0,130.0,135.0,140.0,145.0,150.0 };
            double[] T = new Double[4] { 0.2685,0.7151,0.9534,1.7644 };

            // PutCall identifiers
            string[,] PutCall = new String[NK,NT];
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<=NT-1;t++)
                    PutCall[k,t] = "P";

            // Settings for the Gauthier method
            int time = 1;
            int I1 = 4;
            int I2 = 5;
            double tau = T[time];
            double Put1 = MktPrice[I1,time];
            double Put2 = MktPrice[I2,time];
            double K1 = K[I1];
            double K2 = K[I2];

            // Select arbitrary starting values for kappa, theta, and v0
            double kappa0 = 10;
            double theta0 = 0.1;
            double v00    = 0.01;
            
            // Get Gauthier starting values for rho and sigma
            Gauthier GS = new Gauthier();
            double[] Gstart = GS.GetGauthierValues(kappa0,theta0,v00,settings.S,K1,K2,Put1,Put2,tau,settings.r,settings.q);

            double sigma = Gstart[0];
            double rho   = Gstart[1];

            // Write results
            Console.WriteLine("Gauthier-Possamai starting values");
            Console.WriteLine("---------------------------------");
            Console.WriteLine("Sigma          {0,10:F5}",sigma);
            Console.WriteLine("Rho            {0,10:F5}",rho);
            Console.WriteLine("---------------------------------");
        }
    }
}
