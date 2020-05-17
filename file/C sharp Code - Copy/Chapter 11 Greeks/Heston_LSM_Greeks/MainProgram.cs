using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_LSM_Greeks
{
    class HestonLSMGreeks
    {
        static void Main(string[] args)
        {
            RandomNumbers RN = new RandomNumbers();
            LSMGreeksAlgo LSM = new LSMGreeksAlgo();

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
            HParam param = new HParam();
            param.kappa = 5.0;
            param.theta = 0.16;
            param.sigma = 0.9;
            param.v0 = 0.0625;
            param.rho = 0.1;
            param.lambda = 0.0;

            OpSet settings = new OpSet();
            settings.K = 10.0;
            settings.T = 0.25;
            settings.r = 0.1;
            settings.q = 0.0;
            settings.PutCall = "P";
            settings.trap = 1;

            // Simulation settings
            int NT =  100;
            int NS = 2500;

            // Generate the correlated random variables
            double[,] Zv = new double[NT,NS];
            double[,] Zs = new double[NT,NS];
            double rho = param.rho;
            for(int t=0;t<=NT-1;t++)
            {
                for(int s=0;s<=NS-1;s++)
                {
                    Zv[t,s] = RN.RandomNorm();
                    Zs[t,s] = rho*Zv[t,s] + Math.Sqrt(1-rho*rho)*RN.RandomNorm();
                }
            }

            // Clark and Parrott true prices
            int NK = 5;
            double[] TruePrice = new double[5] { 2.0000,1.107641,0.520030,0.213668,0.082036 };
            double[] S = new double[5] { 8.0, 9.0, 10.0, 11.0, 12.0};

            // Prices and Greeks
            double[] EuroPriceSim = new double[NK];
            double[] EuroDeltaSim = new double[NK];
            double[] EuroGammaSim = new double[NK];
            double[] EuroVega1Sim = new double[NK];
            double[] EuroVannaSim = new double[NK];
            double[] EuroThetaSim = new double[NK];
            double[] EuroRhoSim   = new double[NK];

            double[] AmerPriceSim = new double[NK];
            double[] AmerDeltaSim = new double[NK];
            double[] AmerGammaSim = new double[NK];
            double[] AmerVega1Sim = new double[NK];
            double[] AmerVannaSim = new double[NK];
            double[] AmerThetaSim = new double[NK];
            double[] AmerRhoSim   = new double[NK];

            double[] input = new double[2];

            // Calculate the LSM Euro and Amer prices and Greeks
            for(int k=0;k<=NK-1;k++)
            {
                settings.S = S[k];
                input = LSM.LSMGreeks(settings,param,NT,NS,Zv,Zs,"price");
                EuroPriceSim[k] = input[0];
                AmerPriceSim[k] = input[1];
                input = LSM.LSMGreeks(settings,param,NT,NS,Zv,Zs,"delta");
                EuroDeltaSim[k] = input[0];
                AmerDeltaSim[k] = input[1];
                input = LSM.LSMGreeks(settings,param,NT,NS,Zv,Zs,"gamma");
                EuroGammaSim[k] = input[0];
                AmerGammaSim[k] = input[1];
                input = LSM.LSMGreeks(settings,param,NT,NS,Zv,Zs,"vega1");
                EuroVega1Sim[k] = input[0];
                AmerVega1Sim[k] = input[1];
                input = LSM.LSMGreeks(settings,param,NT,NS,Zv,Zs,"vanna");
                EuroVannaSim[k] = input[0];
                AmerVannaSim[k] = input[1];
                input = LSM.LSMGreeks(settings,param,NT,NS,Zv,Zs,"theta");
                EuroThetaSim[k] = input[0];
                AmerThetaSim[k] = input[1];
                input = LSM.LSMGreeks(settings,param,NT,NS,Zv,Zs,"rho");
                EuroRhoSim[k] = input[0];
                AmerRhoSim[k] = input[1];
            }
            Console.WriteLine("Clarke and Parrott American Greeks with LSM");
            Console.WriteLine("LSM uses {0,3} time steps and {1,3} stock paths",NT,NS);
            Console.WriteLine("------------------------------------------------------------------");
            Console.WriteLine("  S0   Price   Delta    Gamma    Vega1    Vanna    Theta     Rho");
            Console.WriteLine("------------------------------------------------------------------");
            for(int k=0;k<=NK-1;k++)
            {
                Console.WriteLine("{0,3} {1,8:F4} {2,8:F4} {3,8:F4} {4,8:F4} {5,8:F4} {6,8:F4} {7,8:F4}",
                  S[k],AmerPriceSim[k],AmerDeltaSim[k],AmerGammaSim[k],AmerVega1Sim[k],AmerVannaSim[k],AmerThetaSim[k],AmerRhoSim[k]);
            }
            Console.WriteLine("------------------------------------------------------------------");
            Console.WriteLine(" ");
            Console.WriteLine(" ");
            Console.WriteLine("Clarke and Parrott European Greeks with LSM");
            Console.WriteLine("------------------------------------------------------------------");
            Console.WriteLine("  S0   Price   Delta    Gamma    Vega1    Vanna    Theta     Rho");
            Console.WriteLine("------------------------------------------------------------------");
            for(int k=0;k<=NK-1;k++)
            {
                Console.WriteLine("{0,3} {1,8:F4} {2,8:F4} {3,8:F4} {4,8:F4} {5,8:F4} {6,8:F4}",
                  S[k],EuroPriceSim[k],EuroDeltaSim[k],EuroGammaSim[k],EuroVega1Sim[k],EuroVannaSim[k],EuroThetaSim[k],EuroRhoSim[k]);
            }
            Console.WriteLine("----------------------------------------------------------");
        }
    }
}

