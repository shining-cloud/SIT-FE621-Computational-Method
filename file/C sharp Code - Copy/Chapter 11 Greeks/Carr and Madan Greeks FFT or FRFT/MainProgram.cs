using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Carr_and_Madan_FFT_or_FRFT_Greeks
{
    class CarrMadanFFTGreeks
    {
        static void Main(string[] args)
        {
            FastFT FaFT = new FastFT();
            FractionalFFT FrFFT = new FractionalFFT();
            HestonGreeks HG = new HestonGreeks();

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
            param.kappa = 2.0;
            param.theta = 0.06;
            param.sigma = 0.1;
            param.v0 = 0.06;
            param.rho = -0.7;
            param.lambda = 0.0;

            // Option settings
            OpSet settings = new OpSet();
            settings.S = 100.0;
            settings.T = 0.5;
            settings.r = 0.05;
            settings.q = 0.0;
            settings.PutCall = "C";
            settings.trap = 1;

            // FFT settings
            double N2 = 9.0;
            double M = Math.Pow(2.0,N2);
            int N = Convert.ToInt16(M);
            double alpha = 1.5;
            double uplimit = 100.0;
            string rule = "Simpsons";

            // Greeks by the FFT or FRFT
            string method = "FRFT";
            double eta = 0.1;
            double lambdainc = 0.005;
            string[] Greek = { "Price","Delta","Gamma","Theta","Rho","Vega1","Vanna","Volga" };
            double[,] GreeksFFT = new double[N,8];
            double[,] FFT = new double[N,2];
            for(int j=0;j<=7;j++)
            {
                if(method == "FFT")
                    FFT = FaFT.HestonFFTGreek(N,uplimit,alpha,rule,param,settings,Greek[j]);
                else
                    FFT = FrFFT.HestonFRFTGreek(N,uplimit,alpha,rule,lambdainc,eta,param,settings,Greek[j]);
                for(int k=0;k<=N-1;k++)
                    GreeksFFT[k,j] = FFT[k,1];
            }
            double[] K = new double[N];
            for(int k=0;k<=N-1;k++)
                K[k] = FFT[k,0];

            // Select results near the money
            int start = N/2 - 3;
            int end   = N/2 + 3;
            int Nm = end - start + 1;
            double[,] Greeks = new double[Nm,8];
            double[] Strike = new double[Nm];
            for(int k=0;k<=Nm-1;k++)
            {
                Strike[k] = K[start+k];
                for(int j=0;j<=7;j++)
                {
                    Greeks[k,j] = GreeksFFT[start+k,j];
                }
            }

            // Adjust the FFT Greeks vega1 and vanna for v0
            for(int k=0;k<=Nm-1;k++)
            {
                Greeks[k,5] *= 2.0*Math.Sqrt(param.v0);
                Greeks[k,6] *= 2.0*Math.Sqrt(param.v0);
            }

            // Closed-form Greeks near the money and errors
            double[,] GreeksClosed = new double[Nm,8];
            double[] SumAbsError = new double[8];
            for(int j=0;j<=7;j++)
                for (int k=0; k<=Nm-1; k++)
                {
                    settings.K = Strike[k];
                    GreeksClosed[k,j] = HG.CarrMadanGreeks(alpha,settings.T,param,settings,Greek[j],x,w);
                    SumAbsError[j] += Math.Abs(GreeksClosed[k,j] - Greeks[k,j]);
                }

            // Output the results near the money
            Console.WriteLine("-----------------------------------------------------------------------------");
            Console.WriteLine("              Call             Delta              Gamma            Theta");
            Console.WriteLine("Strike   FFT       CM      FFT       CM       FFT       CM     FFT       CM");
            Console.WriteLine("-----------------------------------------------------------------------------");
            for(int k=0;k<=Nm-1;k++)
                Console.WriteLine("{0,6:F2} {1,8:F3} {2,8:F3} {3,8:F3} {4,8:F3} {5,8:F3} {6,8:F3} {7,8:F3} {8,8:F3}",Strike[k],
                    Greeks[k,0],GreeksClosed[k,0],
                    Greeks[k,1],GreeksClosed[k,1],
                    Greeks[k,2],GreeksClosed[k,2],
                    Greeks[k,3],GreeksClosed[k,3]);
            Console.WriteLine("Error {0,15:E3} {1,17:E3} {2,17:E3} {3,17:E3}",SumAbsError[0],SumAbsError[1],SumAbsError[2],SumAbsError[3]);
            Console.WriteLine("-----------------------------------------------------------------------------");
            Console.WriteLine("              Rho              Vega1              Vanna            Volga");
            Console.WriteLine("Strike   FFT       CM      FFT       CM       FFT       CM     FFT       CM");
            Console.WriteLine("-----------------------------------------------------------------------------");
            for(int k=0;k<=Nm-1;k++)
                Console.WriteLine("{0,6:F2} {1,8:F3} {2,8:F3} {3,8:F3} {4,8:F3} {5,8:F3} {6,8:F3} {7,8:F3} {8,8:F3}",Strike[k],
                    Greeks[k,4],GreeksClosed[k,4],
                    Greeks[k,5],GreeksClosed[k,5],
                    Greeks[k,6],GreeksClosed[k,6],
                    Greeks[k,7],GreeksClosed[k,7]);
            Console.WriteLine("Error {0,15:E3} {1,17:E3} {2,17:E3} {3,17:E3}",SumAbsError[4],SumAbsError[5],SumAbsError[6],SumAbsError[7]);
            Console.WriteLine("-----------------------------------------------------------------------------");
        }
    }
}
