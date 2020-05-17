using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_American_Greeks
{
    class MSGreeksHeston
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

            // Settings for Simpson's Rule
            double A = 1.0e-10;
            double B = 150.0;
            int N = 10000;
            int method = 3;

            // Infinite barrier for the European put
            double yinf = 1e4;

            OpSet opset = new OpSet();
            opset.K = Strike;
            opset.r = r;
            opset.q = q;
            opset.T = T;
            opset.PutCall = "P";
            opset.trap = trap;

            // Initialize the barriers, the puts, the pricing errors
            double[] AmerPut = new double[5];
            double[] theta = new double[5];
            double[] y = new double[5];
            double[] error = new double[5];
            double[] Aerror = new double[5];

            // Settings for Bisection algorithm
            double hi = 3.0;
            int MaxIter = 1000;
            double tol = 1.0e-10;
            int NumTerms = 3;

            double[] Price = new double[5];
            double[] Delta = new double[5];
            double[] Gamma = new double[5];
            double[] Vega1 = new double[5];
            double[] Vanna = new double[5];
            double[] Volga = new double[5];
            double[] Theta = new double[5];

            // Find the Medvedev-Scaillet Heston price
            MSGreeks MS = new MSGreeks();
            for(int k=0;k<=4;k++)
            {
                opset.S = S[k];
                Price[k] = MS.MSGreeksFD(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf,"price");
                Delta[k] = MS.MSGreeksFD(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf,"delta");
                Gamma[k] = MS.MSGreeksFD(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf,"gamma");
                Vega1[k] = MS.MSGreeksFD(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf,"vega1");
                Vanna[k] = MS.MSGreeksFD(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf,"vanna");
                Volga[k] = MS.MSGreeksFD(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf,"volga");
                Theta[k] = MS.MSGreeksFD(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf,"theta");
            }

            // Write the results
            Console.WriteLine("--------------------------------------------------------------------------------");
            Console.WriteLine("Medvedev-Scaillet {0:F0}-term approximation to American Greeks",NumTerms);
            Console.WriteLine("Clarke and Parrott prices ");
            Console.WriteLine("--------------------------------------------------------------------------------");
            Console.WriteLine("Spot  TruePrice   Approx     Delta   Gamma    Vega1    Vanna     Volga   Theta");
            Console.WriteLine("--------------------------------------------------------------------------------");
            for(int k=0;k<=4;k++)
                Console.WriteLine("{0,3:F0} {1,10:F5} {2,10:F5} {3,8:F4} {4,8:F4} {5,8:F4} {6,8:F4} {7,8:F4} {8,8:F4}",
                    S[k],TruePrice[k],Price[k],Delta[k],Gamma[k],Vega1[k],Vanna[k],Volga[k],Theta[k]);
            Console.WriteLine("--------------------------------------------------------------------------------");
        }
    }
}
