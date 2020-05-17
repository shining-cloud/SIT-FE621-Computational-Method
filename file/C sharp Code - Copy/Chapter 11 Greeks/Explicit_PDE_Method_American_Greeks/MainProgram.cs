using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Explicit_PDE_Method_American_Greeks
{
    partial class ExplicitPDEAmerican
    {
        static void Main(string[] args)
        {
            // Classes
            Interpolation IP = new Interpolation();
            ExplicitPDE EP = new ExplicitPDE();

            // Number of grid points for the stock, volatility, and maturity
            int nS = 79;       // Stock price
            int nV = 39;        // Volatility
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
            double[] S0 = new double[5] { 8.0,9.0,10.0,11.0,12.0 };
            double[] TruePrice = new double[5] {2.0000, 1.107641, 0.520030, 0.213668, 0.082036};

            // Heston parameters
            HParam param;
            param.kappa =  5;
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
            double Smin = 0.0; double Smax = 2.0*K;
            double Vmin = 0.0; double Vmax = 0.5;
            double Tmin = 0.0; double Tmax = Mat;

            // The maturity time increment and grid
            double dt = (Tmax-Tmin)/Convert.ToDouble(nT);
            double[] T = new double[nT+1];
            for(int i=0;i<=nT;i++)
                T[i] = Convert.ToDouble(i)*dt;

            // Obtain the prices by 2-D interpolation and the errors compared to Clarke and Parrott
            // Also obtain the Greeks
            double V0 = param.v0;
            int N = S0.Length;
            double[] PDEPrice = new double[N];
            double D1,D2,V1,V2,dCdv0,C1,C2,C3,C4,dC2,T1,T2;
            double[] DeltaPDE = new double[N];
            double[] GammaPDE = new double[N];
            double[] Vega1PDE = new double[N];
            double[] VannaPDE = new double[N];
            double[] VolgaPDE = new double[N];
            double[] ThetaPDE = new double[N];
            double dS,dv,c,dz,d,dn;
            double[] z = new double[nS+1];
            double[] S = new double[nS+1];
            double[] n = new double[nV+1];
            double[] V = new double[nV+1];
            for(int k=0;k<=N-1;k++)
            {
                //// Pricing Using a Non-Uniform Grid
                // The stock price grid
                c = S0[k]/5.0;
                dz = 1.0/nS*(IP.aSinh((Smax-S0[k])/c) - IP.aSinh(-S0[k]/c));
                for(int i=0;i<=nS;i++)
                {
                    z[i] = IP.aSinh(-S0[k]/c) + Convert.ToDouble(i)*dz;
                    S[i] = S0[k] + c*Math.Sinh(z[i]);
                }
                S[0] = 0;

                // The volatility grid
                d = Vmax/500.0;
                dn = IP.aSinh(Vmax/d)/nV;
                for(int j=0;j<=nV;j++)
                {
                    n[j] = Convert.ToDouble(j)*dn;
                    V[j] = d*Math.Sinh(n[j]);
                }

                // Solve the PDE and return U(S,v,T) and U(S,v,T-dt);
                Uu output = EP.HestonExplicitPDENonUniformGrid(param,K,r,q,S,V,T,PutCall,EuroAmer);
                double[,] U = output.bigU;
                double[,] u = output.smallU;

                // Price, Delta, Gamma
                dS = 0.01*S0[k];
                dv = 0.01*V0;
                PDEPrice[k] = IP.interp2(V,S,U,V0,S0[k]);
                D1 = IP.interp2(V,S,U,V0,S0[k]+dS);
                D2 = IP.interp2(V,S,U,V0,S0[k]-dS);
                DeltaPDE[k] = (D1-D2)/2.0/dS;
                GammaPDE[k] = (D1-2.0*PDEPrice[k]+D2)/dS/dS;
                // Vega #1
                V1 = IP.interp2(V,S,U,V0+dv,S0[k]);
                V2 = IP.interp2(V,S,U,V0-dv,S0[k]);
                dCdv0 = (V1-V2)/2.0/dv;
                Vega1PDE[k] = dCdv0 * 2.0 * Math.Sqrt(V0);
                // Vanna and volga
                C1 = IP.interp2(V,S,U,V0+dv,S0[k]+dS);
                C2 = IP.interp2(V,S,U,V0-dv,S0[k]+dS);
                C3 = IP.interp2(V,S,U,V0+dv,S0[k]-dS);
                C4 = IP.interp2(V,S,U,V0-dv,S0[k]-dS);
                VannaPDE[k] = (C1-C2-C3+C4)/4.0/dv/dS*2.0*Math.Sqrt(V0);
                dC2 = (V1 - 2*PDEPrice[k] + V2)/dv/dv;
                VolgaPDE[k] = 4.0*Math.Sqrt(V0)*(dC2*Math.Sqrt(V0) + Vega1PDE[k]/4.0/V0);
                // Theta
                T1 = IP.interp2(V,S,U,V0,S0[k]);       // U(S,v,T)
                T2 = IP.interp2(V,S,u,V0,S0[k]);       // U(s,v,T-dt) 
                ThetaPDE[k] = -(T1 - T2)/dt;
            }

            // Output the results
            Console.WriteLine("Stock price grid size  {0}", nS+1);
            Console.WriteLine("Volatility grid size   {0}", nV+1);
            Console.WriteLine("Number of time steps   {0}", nT);
            Console.WriteLine("--------------------------------------------------------------------------");
            Console.WriteLine("Spot    Price      Delta    Gamma    Vega1    Vanna    Volga    Theta");
            Console.WriteLine("--------------------------------------------------------------------------");
            for (int k=0; k<=N-1; k++)
                Console.WriteLine(" {0,2:F0} {1,10:F4} {2,10:F4} {3,8:F4} {4,8:F4} {5,8:F4} {6,8:F4} {7,8:F4}",
                    S0[k],PDEPrice[k],DeltaPDE[k],GammaPDE[k],Vega1PDE[k],VannaPDE[k],VolgaPDE[k],ThetaPDE[k]);
            Console.WriteLine("--------------------------------------------------------------------------");
        }
    }
}


