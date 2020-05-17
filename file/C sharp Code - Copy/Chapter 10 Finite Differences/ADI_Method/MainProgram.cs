using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;
using System.Diagnostics;

namespace ADI_Method
{
    class ADIMethodPricing
    {
        static void Main(string[] args)
        {
            // Classes
            Interpolation IP = new Interpolation();
            HestonPrice HP = new HestonPrice();
            ADIMethod ADI = new ADIMethod();

            // Settings for the clock;
            Stopwatch sw = new Stopwatch();
            TimeSpan ts = sw.Elapsed;

            // Settings for the option price calculation
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

            // Strike price, risk free rate, dividend yield, and maturity
            double K   = 100.0;
            double r   = 0.02;
            double q   = 0.05;
            double Mat = 0.15;

            // Heston parameters.  Case 1 of Hout and Foulon (Table 1)
            HParam param;
            param.kappa =  1.50;
            param.theta =  0.04;
            param.sigma =  0.30;
            param.rho   = -0.90;
            param.v0    =  0.05;
            param.lambda = 0.00;

            // Minimum and maximum values for the Stock Price, Volatility, and Maturity
            double Smin = 0.0; double Smax = 2.0*K;
            double Vmin = 0.0; double Vmax = 0.5;
            double Tmin = 0.0; double Tmax = Mat;

            double[] S,V;

            // Points for stock, vol, maturity
            int nS = 9;
            int nV = 9;
            int nT = 9;
            int NS,NV,NT;

            // Select the grid type
            string GridType = "NonUniform";

            if(GridType == "Uniform")
            {
                // Increment for Stock Price, Volatility, and Maturity
                double ds = (Smax-Smin)/Convert.ToDouble(nS);
                double dv = (Vmax-Vmin)/Convert.ToDouble(nV);

                // Grid Vectors for the Stock Price, Volatility, and Maturity
                NS = nS+1;
                NV = nV+1;
                S = new double[NS];
                V = new double[NV];
                for(int s=0;s<=NS-1;s++)
                    S[s] = Convert.ToDouble(s)*ds;
                for(int v=0;v<=NV-1;v++)
                    V[v] = Convert.ToDouble(v)*dv;
            }
            else
            {
                // The stock price grid
                double c = K/5.0;
                double dz = 1.0/nS*(IP.aSinh((Smax-K)/c) - IP.aSinh(-K/c));
                double[] z = new double[nS+1];
                S = new double[nS+1];
                for(int i=0;i<=nS;i++)
                {
                    z[i] = IP.aSinh(-K/c) + Convert.ToDouble(i)*dz;
                    S[i] = K + c*Math.Sinh(z[i]);
                }
                S[0] = 0.0;

                // The volatility grid
                double d = Vmax/10.0;
                double dn = IP.aSinh(Vmax/d)/nV;
                double[] n = new double[nV+1];
                V = new double[nV+1];
                for(int j=0;j<=nV;j++)
                {
                    n[j] = Convert.ToDouble(j)*dn;
                    V[j] = d*Math.Sinh(n[j]);
                }
                NS = nS+1;
                NV = nV+1;
            }

            // The maturity time increment and grid
            NT = nT+1;
            double[] T = new double[NT];
            double dt = (Tmax-Tmin)/Convert.ToDouble(nT);
            for(int i=0;i<=nT;i++)
                T[i] = Convert.ToDouble(i)*dt;

            // Settings for the option price
            double thet;
            double S0 = 101.52;
            double V0 = 0.05412;

            // Obtain the Exact Price
            string PutCall = "C";
            int trap = 1;
            param.v0 = V0;
            double HPrice = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);

            // Obtain the ADI Prices and errors
            double[] EPrice = new double[4];
            double[] IPrice = new double[4];
            double[] CPrice = new double[4];
            double[] EError = new double[4];
            double[] IError = new double[4];
            double[] CError = new double[4];

            string[] scheme = new string[4] { "DO","CS","MCS","HV" };

            for(int s=0;s<=3;s++)
            {
                sw.Reset();
                sw.Start();

                // Obtain the ADI price with theta = 0 (Explicit)
                thet = 0.0;
                EPrice[s] = ADI.ADIPrice(scheme[s],thet,param,S0,V0,K,r,q,S,V,T,GridType);
                EError[s] = EPrice[s] - HPrice;

                // Obtain the ADI price with theta = 1 (Implicit)
                thet = 1.0;
                IPrice[s] = ADI.ADIPrice(scheme[s],thet,param,S0,V0,K,r,q,S,V,T,GridType);
                IError[s] = IPrice[s] - HPrice;

                // Obtain the ADI price with theta = 0.5 (Crank Nicolson)
                thet = 0.5;
                CPrice[s] = ADI.ADIPrice(scheme[s],thet,param,S0,V0,K,r,q,S,V,T,GridType);
                CError[s] = CPrice[s] - HPrice;

                sw.Stop();
                ts = sw.Elapsed;
                Console.WriteLine("The {0,3:0} scheme ADI prices were calculated in {1:0} minutes and {2,5:F3} seconds ",scheme[s],ts.Minutes,ts.Seconds);
            }

            // Output the results
            Console.WriteLine("-----------------------------------------------------------------------");
            Console.WriteLine("Stock price grid size of {0:0} ",NS);
            Console.WriteLine("Volatility grid size of  {0:0} ",NV);
            Console.WriteLine("Number of time steps     {0:0} ",NT);
            Console.WriteLine("Closed form Price is     {0,5:F6}",HPrice);
            Console.WriteLine("-----------------------------------------------------------------------");
            Console.WriteLine("ADI    Explicit   Error    Implicit   Error  Crank-Nicolson  Error");
            Console.WriteLine("-----------------------------------------------------------------------");
            for(int s=0;s<=3;s++)
                Console.WriteLine("{0,3} {1,10:F4} {2,8:F4} {3,10:F4} {4,8:F4} {5,10:F4} {6,11:F4}",scheme[s],EPrice[s],EError[s],IPrice[s],IError[s],CPrice[s],CError[s]);
            Console.WriteLine("-----------------------------------------------------------------------");
        }
    }
}
