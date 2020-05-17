using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Explicit_PDE_Method
{
    class ExplicitPDE
    {
        static void Main(string[] args)
        {
            // Classes
            PDEAlgoNU NUG = new PDEAlgoNU();
            PDEAlgoU UG = new PDEAlgoU();
            Interpolation IP = new Interpolation();
            HestonPrice HP = new HestonPrice();

            // Illustration of pricing using uniform and non-uniform grids
            // Strike price, risk free rate, dividend yield, and maturity
            double K = 100.0;
            double r = 0.02;
            double q = 0.05;
            double Mat = 0.15;

            // Heston parameters
            HParam param;
            param.kappa =  1.5;
            param.theta =  0.04;
            param.sigma =  0.3;
            param.rho   = -0.9;
            param.v0    =  0.05;
            param.lambda = 0;

            // Minimum and maximum values for the Stock Price, Volatility, and Maturity
            double Smin = 0.0; double Smax = 2.0*K;
            double Vmin = 0.0; double Vmax = 0.5;
            double Tmin = 0.0; double Tmax = Mat;

            // Number of grid points for the stock, volatility, and maturity
            int nS = 79;        // Stock price
            int nV = 39;        // Volatility
            int nT = 3000;      // Maturity

            // The maturity time increment and grid
            double dt = (Tmax-Tmin)/Convert.ToDouble(nT);
            double[] T = new double[nT+1];
            for(int i=0;i<=nT;i++)
                T[i] = Convert.ToDouble(i)*dt;

            //// Pricing Using a Uniform Grid
            // Increment for Stock Price and volatility
            double ds = (Smax-Smin)/Convert.ToDouble(nS);
            double dv = (Vmax-Vmin)/Convert.ToDouble(nV);

            // The stock price and volatility grids
            double[] S = new double[nS+1];
            for(int i=0;i<=nS;i++)
                S[i] = Convert.ToDouble(i)*ds;

            double[] V = new double[nV+1];
            for(int i=0;i<=nV;i++)
                V[i] = Convert.ToDouble(i)*dv;

            // Solve the PDE
            double[,] U = UG.HestonExplicitPDE(param,K,r,q,S,V,T);

            // Obtain the price by 2-D interpolation
            double S0 = 101.52;
            double V0 = 0.05412;
            double UniformPrice = IP.interp2(V,S,U,V0,S0);

            // Pricing Using a Non-Uniform Grid
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
            S[0] = 0;

            // The volatility grid
            double d = Vmax/500.0;
            double dn = IP.aSinh(Vmax/d)/nV;
            double[] n = new double[nV+1];
            V = new double[nV+1];
            for(int j=0;j<=nV;j++)
            {
                n[j] = Convert.ToDouble(j)*dn;
                V[j] = d*Math.Sinh(n[j]);
            }

            // Solve the PDE
            double[,] UU = NUG.HestonExplicitPDENonUniformGrid(param,K,r,q,S,V,T);

            // Obtain the price by 2-D interpolation
            double NonUniformPrice = IP.interp2(V,S,UU,V0,S0);

            //// True Price
            int trap = 1;
            string PutCall = "C";

            // Values at which to obtain the true price
            param.v0 = V0;

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

            // The closed form price
            double ClosedPrice = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);

            // The errors
            double UError = UniformPrice - ClosedPrice;
            double NError = NonUniformPrice - ClosedPrice;

            // Output the results
            Console.WriteLine("----------------------------------------------");
            Console.WriteLine("Grid sizes");
            Console.WriteLine("  Stock price: {0:0}, Volatility: {1:0}, Time: {2:0}  ",nS+1,nV+1,nT);
            Console.WriteLine("----------------------------------------------");
            Console.WriteLine("  Method              Price        Error");
            Console.WriteLine("----------------------------------------------");
            Console.WriteLine("  Closed Form        {0,5:F4} ", ClosedPrice);
            Console.WriteLine("  Uniform Grid       {0,5:F4}    {1,10:F4}",UniformPrice,UError);
            Console.WriteLine("  Non-Uniform Grid   {0,5:F4}    {1,10:F4}",NonUniformPrice,NError);
            Console.WriteLine("----------------------------------------------");
            Console.WriteLine();
        }
    }
}



