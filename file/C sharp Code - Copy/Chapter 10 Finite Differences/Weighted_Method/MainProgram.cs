using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;
using System.Diagnostics;

namespace Weighted_Method
{
    class WeightedMethod
    {
        static void Main(string[] args)
        {
            // Classes
            Interpolation IP = new Interpolation();
            BuildDerivativesNU BNU = new BuildDerivativesNU();
            BuildDerivativesU BU = new BuildDerivativesU();
            HestonPrice HP = new HestonPrice();
            MatrixOps MO = new MatrixOps();
            WeightedPriceAlgo WP = new WeightedPriceAlgo();

            // Settings for the option price calculation
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] X = new double[32];
            double[] W = new double[32];
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

            // Settings for the option price
            double S0 = 101.52;
            double V0 = 0.05412;

            // Heston parameters.  Case 1 of Hout and Foulon (Table 1)
            HParam param;
            param.kappa =  1.50;
            param.theta =  0.04;
            param.sigma =  0.30;
            param.rho   = -0.90;
            param.v0    =  V0;
            param.lambda = 0.00;

            // Obtain the Exact Price
            string PutCall = "C";
            int trap = 1;
            double HPrice = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);

            // Minimum and maximum values for the Stock Price, Volatility, and Maturity
            double Smin = 0.0; double Smax = 2.5*K;
            double Vmin = 0.0; double Vmax = 0.5;
            double Tmin = 0.0; double Tmax = Mat;

            // Points for stock, vol, maturity
            double[] S,V;
            int nS = 29;
            int nV = 29;
            int nT = 19;
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
            else //if (GridType == "NonUniform")
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


            // Obtain the submatrices of L
            LMatrices LMat;
            if (GridType == "Uniform")
                LMat = BU.BuildDerivatives(S,V,T);
            else
                LMat = BNU.BuildDerivativesNonUniform(S,V,T);
            double[,] derS = LMat.derS;
            double[,] derSS = LMat.derSS;
            double[,] derV1 = LMat.derV1;
            double[,] derV2 = LMat.derV2;
            double[,] derVV = LMat.derVV;
            double[,] derSV = LMat.derSV;
            double[,] R     = LMat.R;

            // Build the L matrix
            int N = NS*NV;
            double[,] L = new double[N,N];
            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=N-1;j++)
                    L[i,j] = (r-q)*derS[i,j] + param.kappa*param.theta*derV1[i,j] - param.kappa*derV2[i,j] + 0.5*derSS[i,j] 
                           + 0.5*param.sigma*param.sigma*derVV[i,j] + param.rho*param.sigma*derSV[i,j] - r*R[i,j];

            // Identity and A and B matrix
            double[,] I = new double[N,N];
            double[,] A = new double[N,N];
            double[,] B = new double[N,N];
            double[,] invA = new double[N,N];
            for(int i=0;i<=N-1;i++)
                I[i,i] = 1.0;

            // Settings for the clock
            Stopwatch sw = new Stopwatch();
            TimeSpan ts = sw.Elapsed;

            // Obtain the Explicit, Implicit, and Crank-Nicolson prices
            double[] thet = { 0.0,1.0,0.5 };
            double[] Price = new double[3];
            double[] Error = new double[3];
            for(int k=0;k<=2;k++)
            {
                for(int i=0;i<=N-1;i++)
                    for(int j=0;j<=N-1;j++)
                    {
                        A[i,j] = I[i,j] - thet[k]*dt*L[i,j];
                        B[i,j] = I[i,j] + (1.0-thet[k])*dt*L[i,j];
                    }
                sw.Reset();
                sw.Start();
                if(thet[k] == 0.0)
                    invA = MO.CreateI(N);
                else
                    invA = MO.MInvLU(A);
                ts = sw.Elapsed;
                Console.WriteLine("Calculated the inverse in {0:0}-{1:0}-{2:0} min-sec-msec",ts.Minutes,ts.Seconds,ts.Milliseconds);
                Price[k] = WP.WeightedPrice(thet[k],L,S0,V0,K,r,q,Mat,S,V,T,A,invA,B);
                Error[k] = Price[k] - HPrice;
            }
            // Output the results
            Console.WriteLine("-------------------------------------------------");
            Console.WriteLine("Grid type    : {0}",GridType);
            Console.WriteLine("-------------------------------------------------");
            Console.WriteLine("Stock price grid size of {0:0} ",NS);
            Console.WriteLine("Volatility grid size of  {0:0} ",NV);
            Console.WriteLine("Number of time steps     {0:0} ",NT);
            Console.WriteLine("-------------------------------------------------");
            Console.WriteLine("Method                    Price      Dollar Error");
            Console.WriteLine("-------------------------------------------------");
            Console.WriteLine("Closed form              {0,5:F4}",HPrice);
            Console.WriteLine("Explicit method          {0,5:F4} {1,10:F4}",Price[0],Error[0]);
            Console.WriteLine("Implicit method          {0,5:F4} {1,10:F4}",Price[1],Error[1]);
            Console.WriteLine("Crank Nicolson           {0,5:F4} {1,10:F4}",Price[2],Error[2]);
            Console.WriteLine("-------------------------------------------------");
        }
    }
}
