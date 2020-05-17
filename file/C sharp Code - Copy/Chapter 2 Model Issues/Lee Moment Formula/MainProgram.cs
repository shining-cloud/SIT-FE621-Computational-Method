using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

// Heston parameters
public struct HParam
{
    public double kappa;        // Mean reversion speed
    public double theta;        // Mean reversion level
    public double sigma;        // Volatility of variance
    public double v0;           // Initial variance
    public double rho;          // Correlation
    public double lambda;       // Risk parameter
}

namespace Lee_Moment_Formula
{
    class LeeMomentCalculations
    {
        static void Main(string[] args)
        {
            // Underlying settings
            double S = 129.14;
            double tau = 37.0/365.0;
            double r = 0.0;
            double q = 0.0;

            // Parameter settings
            HParam param;
            param.kappa =  4.6987;
            param.theta =  0.1250;
            param.sigma =  2.1149;
            param.v0    =  0.0200;
            param.rho   = -0.2741;
            param.lambda = 0.0;
            int trap = 1;

            // Moment explosion setting
            int HiLimit =  8;
            int LoLimit = -5;

            // Find the slopes and parameter bounds
            LeeMoment LM = new LeeMoment();
            double[] output = LM.FindLeeBounds(S,r,q,tau,param,trap,LoLimit,HiLimit);
            double bR = output[0];
            double bL = output[1];
            double LowerAP = output[2];
            double UpperAP = output[3];

            // Output the results
            Console.WriteLine("Andersen and Piterbarg exploding moment bounds");
            Console.WriteLine("----------------------------------------------");
            Console.WriteLine("Upper Bound {0,10:F5}",UpperAP);
            Console.WriteLine("Lower Bound {0,10:F5}",LowerAP);
            Console.WriteLine(" ");
            Console.WriteLine("Roger Lee Moment formulas");
            Console.WriteLine("----------------------------------------------");
            Console.WriteLine("Right Tail Slope {0,10:F5}",bR);
            Console.WriteLine("Left Tail Slope  {0,10:F5}",bL);
            Console.WriteLine(" ");
        }
    }
}


