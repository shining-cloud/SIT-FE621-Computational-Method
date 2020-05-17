using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Explicit_PDE_Method_American
{
    class ExplicitPDEPrice
    {
        static void Main(string[] args)
        {
            // Illustration of American option pricing using non-uniform grids

            // Number of grid points for the stock, volatility, and maturity
            int nS = 139;        // Stock price
            int nV = 89;        // Volatility
            int nT = 3000;      // Maturity

            // Option flavor
            string PutCall = "P";
            string EuroAmer = "A";

            // Strike price, risk free rate, dividend yield, and maturity
            // True prices from Clarke and Parrott (1999)
            double K = 10.0;
            double r = 0.1;
            double q = 0.00;
            double Mat = 0.25;
            double[] Spot = new double[5] { 8.0,9.0,10.0,11.0,12.0 };
            double[] TruePrice = new double[5] {2.0, 1.107641, 0.520030, 0.213668, 0.082036};

            // Heston parameters
            HParam param;
            param.kappa =  5.0;
            param.theta =  0.16;
            param.sigma =  0.9;
            param.rho   =  0.1;
            param.v0    =  0.0625;
            param.lambda = 0;

            // Settings for the European option price calculation
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

            // Minimum and maximum values for the Stock Price, Volatility, and Maturity
            double Smin = 0.0; double Smax = 3.0*K;
            double Vmin = 0.0; double Vmax = 0.5;
            double Tmin = 0.0; double Tmax = Mat;

            // The maturity time increment and grid
            double dt = (Tmax-Tmin)/Convert.ToDouble(nT);
            double[] T = new double[nT+1];
            for(int i=0;i<=nT;i++)
                T[i] = Convert.ToDouble(i)*dt;

            // Pricing Using a Non-Uniform Grid
            // The stock price grid
            Interpolation IP = new Interpolation();
            double c = K/1.0;
            double dz = 1.0/nS*(IP.aSinh((Smax-K)/c) - IP.aSinh(-K/c));
            double[] z = new double[nS+1];
            double[] S = new double[nS+1];
            for(int i=0;i<=nS;i++)
            {
                z[i] = IP.aSinh(-K/c) + Convert.ToDouble(i)*dz;
                S[i] = K + c*Math.Sinh(z[i]);
            }
            S[0] = 0;

            // The volatility grid
            double d = Vmax/10.0;
            double dn = IP.aSinh(Vmax/d)/nV;
            double[] n = new double[nV+1];
            double[] V = new double[nV+1];
            for(int j=0;j<=nV;j++)
            {
                n[j] = Convert.ToDouble(j)*dn;
                V[j] = d*Math.Sinh(n[j]);
            }

            // Solve the PDE
            ExplicitPDE EP = new ExplicitPDE();
            double[,] UU = EP.HestonExplicitPDENonUniformGrid(param,K,r,q,S,V,T,PutCall,EuroAmer);

            // Obtain the prices by 2-D interpolation and the errors compared to Clarke and Parrott
            // Obtain also the European price
            Hestonprice HP = new Hestonprice();
            int trap = 1;
            double V0 = param.v0;
            double S0;
            int N = Spot.Length;
            double[] AmerPrice = new double[N];
            double[] Error = new double[N];
            double[] EuroPrice = new double[N];
            for(int k=0;k<=N-1;k++)
            {
                S0 = Spot[k];
                AmerPrice[k] = IP.interp2(V,S,UU,V0,S0);
                Error[k] = TruePrice[k] - AmerPrice[k];
                EuroPrice[k] = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);
            }

            // Output the results
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine("Explicit PDE method for American puts - Clarke and Parrot (1999) example");
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine("  Grid sizes");
            Console.WriteLine("  Stock price: {0:0}, Volatility: {1:0}, Time: {2:0}  ",nS+1,nV+1,nT);
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine("  Strike EuroPrice  TruePrice   PDEAmerPrice     Error");
            Console.WriteLine("------------------------------------------------------------------------");
            for(int k=0;k<=N-1;k++)
                Console.WriteLine("  {0,3:0}    {1,8:F6}    {2,8:F6}    {3,8:F6}  {4,12:F6}",Spot[k],EuroPrice[k],TruePrice[k],AmerPrice[k],Error[k]);
            Console.WriteLine("------------------------------------------------------------------------");
            Console.WriteLine();
        }
    }
}


