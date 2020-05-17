using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

// Heston parameters
public struct HParam
{
    public double kappa;        // Mean reversion speed
    public double theta;        // Mean reversion level
    public double sigma;        // Volatility of variance
    public double v0;           // Initial variance
    public double rho;          // Correlation
}

namespace Variance_Swap
{
    class VarianceSwap
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] x = new Double[32];
            double[] w = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    x[k] = double.Parse(bits[0]);
                    w[k] = double.Parse(bits[1]);
                }

            // Strikes
            double[] K = new double[26] {1100.0,1105.0,1110.0,1115.0,1120.0,1125.0,1130.0,
                                         1135.0,1140.0,1145.0,1150.0,1155.0,1160.0,1165.0,
                                         1170.0,1175.0,1180.0,1185.0,1190.0,1195.0,1200.0,
                                         1205.0,1210.0,1215.0,1220.0,1225.0};
            int NK = K.Length;

            // Settings for the underlying
            double S  = 1164.97;
            double r = 0.0;
            double q  = 0.0071;
            double T = 15.0/365.0;

            // Parameter estimates
            HParam param;
            param.kappa =  3.0002;
            param.theta =  0.0633;
            param.sigma =  0.5015;
            param.v0    =  0.1044;
            param.rho   = -0.9008;

            // Heston variance swap
            double Hvar = (param.v0-param.theta)*(1.0-Math.Exp(-param.kappa*T))/param.kappa/T + param.theta;

            // Demeterfi et al variance swap fair strike
            // Implied volatility settings
            double a = 0.001;
            double b = 10.0;
            double Tol = 1e-10;
            int MaxIter = 1000;
            int trap = 1;

            // Generate the calls and puts and find their implied volatilities
            double[] Call = new double[NK];
            double[] Put  = new double[NK];
            double[] CallIV = new double[NK];
            double[] PutIV  = new double[NK];

            HestonPrice HP = new HestonPrice();
            Bisection B = new Bisection();
            for(int k=0;k<=NK-1;k++)
            {
                Call[k] = HP.HestonPriceGaussLaguerre(param,S,K[k],r,q,T,trap,"C",x,w);
                Put[k] = HP.HestonPriceGaussLaguerre(param,S,K[k],r,q,T,trap,"P",x,w);
                CallIV[k] = B.BisecBSIV("C",S,K[k],r,q,T,a,b,Call[k],Tol,MaxIter);
                 PutIV[k] = B.BisecBSIV("P",S,K[k],r,q,T,a,b,Put[k],Tol,MaxIter);
            }

            // Replication algorithm for the variance swap
            // Index for ATM
            int ATM = 13;

            // Number of call prices to create by interpolation
            int NI = 2000;

            // OTM calls and their IV
            double[] KC = new double[NK-ATM];
            double[] CallV = new double[NK-ATM];
            for(int k=0;k<=NK-ATM-1;k++)
            {
                KC[k] = K[ATM+k];
                CallV[k] = CallIV[ATM+k];
            }
            // Interpolation of strikes and OTM calls
            Interpolate I = new Interpolate();
            double dKC = (KC[NK-ATM-1] - KC[0])/Convert.ToDouble(NI);
            double[] KCI = new double[NI+1];
            double[] CallVI = new double[NI+1];
            for(int k=0;k<=NI;k++)
            {
                KCI[k] = KC[0] + dKC*Convert.ToDouble(k);
                CallVI[k] = I.interp1(KC,CallV,KCI[k]);
            }

            // OTM puts and their IV
            double[] KP = new double[ATM+1];
            double[] PutV = new double[ATM+1];
            for(int k=0;k<=ATM;k++)
            {
                KP[k] = K[k];
                PutV[k] = PutIV[k];
            }
            // Interpolation of strikes and OTM puts
            double dKP = (KP[ATM] - KP[0])/Convert.ToDouble(NI);
            double[] KPI = new double[NI+1];
            double[] PutVI = new double[NI+1];
            for(int k=0;k<=NI;k++)
            {
                KPI[k] = KP[0] + dKP*Convert.ToDouble(k);
                PutVI[k] = I.interp1(KP,PutV,KPI[k]);
            }

            // Calculate the Demeterfi et al (1999) fair strike
            VarSwap VS = new VarSwap();
            double Kvar = VS.VarianceSwap(KCI,CallVI,KPI,PutVI,S,T,r,q);

            // Output the results
            Console.WriteLine("Derman et al variance swap fair strike");
            Console.WriteLine("Using {0} Interpolation points",NI);
            Console.WriteLine("---------------------------------------------------------------");
            Console.WriteLine("Estimate of fair volatility using replication       {0:F8}",Kvar);
            Console.WriteLine("Estimate of fair volatility using Heston parameters {0:F8}",Hvar);
            Console.WriteLine("---------------------------------------------------------------");
        }
    }
}
