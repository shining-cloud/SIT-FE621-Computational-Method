using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Weighted_Method
{
    class WeightedPriceAlgo
    {
        public double WeightedPrice(double thet,double[,] L,double S0,double V0,double K,double r,double q,double Mat,double[] S,double[] V,double[] T,double[,] A,double[,] invA,double[,] B)
        {
            // Heston Call price using the Weighted Method
            // Requires a uniform grid for the stock price, volatility, and maturity
            // INPUTS
            //   thet = theta parameter for the Weighted method
            //   L    = operator matrix for the Heston model
            //   params = vector of Heston parameters
            //   S0 = Spot price at which to price the call
            //   V0 = variance at which to price the call
            //   K  = Strike price
            //   r  = risk free rate
            //   q  = dividend yield
            //   Mat = maturity
            //   S = uniform grid for the stock price
            //   V = uniform grid for the volatility
            //   T = uniform grid for the maturity
            //   A = "A" matrix for weighted method
            //   invA = A inverse
            //   B = "B" matrix for weighted method
            // OUTPUT
            //   y = 2-D interpolated European Call price

            MatrixOps MO = new MatrixOps();
            Interpolation IP = new Interpolation();

            // Required vector lengths and time increment
            int NS = S.Length;
            int NV = V.Length;
            int NT = T.Length;
            double dt = T[1]-T[0];

            // Identity matrix
            int N = NS*NV;
            double[,] I = MO.CreateI(N);

            // Initialize the U and u vectors
            double[] U = new double[N];
            double[] u = new double[N];

            // U(0) vector - value of U(T) at maturity
            double[] Si = new double[N];
            int k = 0;
            for(int v=0;v<=NV-1;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    Si[k] = S[s];
                    U[k] = Math.Max(Si[k] - K, 0.0);
                    k += 1;
                }

            // Loop through the time increments, updating U(t) to U(t+1) at each step
            for (int t=2; t<=NT; t++)
            {
                for (k=0; k<=N-1; k++)
                    u[k] = U[k];
                if(thet==0.0)
                    U = MO.MVMult(B,u);                    // Explicit Method
                else if(thet==1.0)
                    U = MO.MVMult(invA,u);                 // Implicit Method
                else
                    U = MO.MVMult(invA,MO.MVMult(B,u));       // Weighted Methods
            }

            // Restack the U vector to output a matrix
            double[,] UU = new double[NS,NV];
            k = 0;
            for(int v=0;v<=NV-1;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    UU[s,v] = U[k];
                    k += 1;
                }

            // Interpolate to get the price at S0 and v0
            return IP.interp2(V,S,UU,V0,S0);
        }
    }
}



