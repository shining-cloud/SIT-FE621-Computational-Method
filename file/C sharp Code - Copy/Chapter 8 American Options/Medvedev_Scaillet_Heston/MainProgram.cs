using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Heston
{
    partial class MS
    {
        static void Main(string[] args)
        {
            // Settings from Clarke and Parrott
            double[] S = new double[5] { 8.0,9.0,10.0,11.0,12.0 };
            double[] TruePrice = new double[5] { 2.00,1.107641,0.520030,0.213668,0.082036 };
            double Strike = 10;
            double r = 0.1;
            double q = 0;
            double T = 0.25;

            // Heston parameters
            double kappaV = 5.0;
            double thetaV = 0.16;
            double sigmaV = 0.9;
            double rho = 0.1;
            double v0  = 0.0625;
            double lambda = 0.0;
            int trap = 1;
            HParam param = new HParam();
            param.kappa = kappaV;
            param.theta = thetaV;
            param.sigma = sigmaV;
            param.v0 = v0;
            param.rho = rho;
            param.lambda = lambda;

            // Option settings
            OpSet opset = new OpSet();
            opset.K = Strike;
            opset.r = r;
            opset.q = q;
            opset.T = T;
            opset.PutCall = "P";
            opset.trap = trap;

            // Settings for Medvedev-Scaillet expansion
            MSSet msset = new MSSet();
            msset.yinf = 1e4;
            msset.A = 1.0e-10;
            msset.B = 150.0;
            msset.N = 10000;
            msset.method = 3;
            msset.a = -25.0;
            msset.b = +25.0;
            msset.MaxIter = 50000;
            msset.tol = 1.0e-10;
            msset.dt = 1e-5;
            msset.NumTerms = 4;

            // Initialize the barriers, the puts, the pricing errors
            double[] MSAmerPut = new double[5];
            double[] theta = new double[5];
            double[] y = new double[5];
            double[] error = new double[5];
            double[] Aerror = new double[5];
            double[] output = new double[6];
            double[] Theta = new double[5];

            // Find the Medvedev-Scaillet Heston price
            MSPriceHeston MS = new MSPriceHeston();
            for(int k=0;k<=4;k++)
            {
                opset.S = S[k];
                output = MS.MSPrice(param,opset,msset);
                MSAmerPut[k] = output[2];
                theta[k]   = output[4];
                y[k]       = output[5];
                error[k] = TruePrice[k] - MSAmerPut[k];
                Aerror[k] = Math.Abs(error[k]);
            }
            double TotalError = Aerror.Sum();

            // Write the results
            Console.WriteLine("-----------------------------------------------------------");
            Console.WriteLine("Medvedev-Scaillet {0:F0}-term approximation",msset.NumTerms);
            Console.WriteLine("Clarke and Parrott prices ");
            Console.WriteLine("-----------------------------------------------------------");
            Console.WriteLine("Spot  TruePrice   MSPrice    Error      Barrier   Moneyness");
            Console.WriteLine("-----------------------------------------------------------");
            for(int k=0;k<=4;k++)
                Console.WriteLine("{0,3:F0} {1,10:F5} {2,10:F5} {3,10:F5} {4,10:F5} {5,10:F5}",S[k],TruePrice[k],MSAmerPut[k],error[k],y[k],theta[k]);
            Console.WriteLine("-----------------------------------------------------------");
            Console.WriteLine("Total Absolute Error with {0:F0} terms {1:F5}",msset.NumTerms,TotalError);
            Console.WriteLine("-----------------------------------------------------------");
        }
    }
}

