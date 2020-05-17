using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Explicit_PDE_Method_Greeks
{
    partial class ExplicitPDEGreeks
    {
        static void Main(string[] args)
        {
            // Classes
            Interpolation IP = new Interpolation();
            ExplicitPDE EP = new ExplicitPDE();
            HestonPrice HP = new HestonPrice();

            // Illustration of pricing using uniform and non-uniform grids
            // Strike price, risk free rate, dividend yield, and maturity
            double K = 100.0;
            double r = 0.02;
            double q = 0.0;
            double Mat = 0.15;

            // Heston parameters
            HParam param;
            param.kappa =  1.5;
            param.theta =  0.04;
            param.sigma =  0.3;
            param.rho   = -0.9;
            param.v0    =  0.05;
            param.lambda = 0.0;

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

            //// Pricing Using a Non-Uniform Grid
            // The stock price grid
            double c = K/5.0;
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
            double d = Vmax/500.0;
            double dn = IP.aSinh(Vmax/d)/nV;
            double[] n = new double[nV+1];
            double[] V = new double[nV+1];
            for(int j=0;j<=nV;j++)
            {
                n[j] = Convert.ToDouble(j)*dn;
                V[j] = d*Math.Sinh(n[j]);
            }

            // Solve the PDE and return U(S,v,T) and U(S,v,T-dt)
            Uu output = EP.HestonExplicitPDENonUniformGrid(param,K,r,q,S,V,T);
            double[,] U = output.bigU;
            double[,] u = output.smallU;

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

            // Values
            double S0 = 101.52;
            double V0 = 0.05412;
            int trap = 1;
            string PutCall = "C";

            // The closed form price and finite difference Greeks
            param.v0 = V0;
            double PriceClosed = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);
            double dS = 1.0;
            double dV = 1e-2;

            // Delta
            double D1 = HP.HestonPriceGaussLaguerre(param,S0+dS,K,r,q,Mat,trap,PutCall,X,W);
            double D2 = HP.HestonPriceGaussLaguerre(param,S0-dS,K,r,q,Mat,trap,PutCall,X,W);
            double DeltaFD = (D1-D2)/2.0/dS;

            // Gamma
            double GammaFD = (D1 - 2.0*PriceClosed + D2)/dS/dS;

            // Vega #1
            param.v0 = V0+dV;
            double V1 = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);
            param.v0 = V0-dV;
            double V2 = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);
            double Vega1FD = (V1-V2)/2.0/dV*2.0*Math.Sqrt(V0);

            // Vanna
            param.v0 = V0+dV;
            double C1 = HP.HestonPriceGaussLaguerre(param,S0+dS,K,r,q,Mat,trap,PutCall,X,W);
            double C3 = HP.HestonPriceGaussLaguerre(param,S0-dS,K,r,q,Mat,trap,PutCall,X,W);
            param.v0 = V0-dV;
            double C2 = HP.HestonPriceGaussLaguerre(param,S0+dS,K,r,q,Mat,trap,PutCall,X,W);
            double C4 = HP.HestonPriceGaussLaguerre(param,S0-dS,K,r,q,Mat,trap,PutCall,X,W);
            double VannaFD = (C1 - C2 - C3 + C4)/4.0/dV/dS*2.0*Math.Sqrt(V0);
            param.v0 = V0;
	
            // Volga
            double dC2 = (V1 - 2.0*PriceClosed + V2) /(dV*dV);
            double VolgaFD = 4.0*Math.Sqrt(V0)*(dC2*Math.Sqrt(V0) + Vega1FD/4.0/V0);

            // Theta
            double T1 = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat,trap,PutCall,X,W);
            double T2 = HP.HestonPriceGaussLaguerre(param,S0,K,r,q,Mat-dt,trap,PutCall,X,W);
            double ThetaFD = -(T1-T2)/dt;

            // The PDE price and Greeks
            double PricePDE = IP.interp2(V,S,U,V0,S0);

            // Delta
            D1 = IP.interp2(V,S,U,V0,S0+dS);
            D2 = IP.interp2(V,S,U,V0,S0-dS);
            double DeltaPDE = (D1-D2)/2.0/dS;

            // Gamma
            double GammaPDE = (D1 - 2.0*PricePDE + D2)/dS/dS;

            // Vega #1
            V1 = IP.interp2(V,S,U,V0+dV,S0);
            V2 = IP.interp2(V,S,U,V0-dV,S0);
            double Vega1PDE = (V1-V2)/2.0/dV*2.0*Math.Sqrt(V0);

            // Vanna
            C1 = IP.interp2(V,S,U,V0+dV,S0+dS);
            C2 = IP.interp2(V,S,U,V0-dV,S0+dS);
            C3 = IP.interp2(V,S,U,V0+dV,S0-dS);
            C4 = IP.interp2(V,S,U,V0-dV,S0-dS);
            double VannaPDE = (C1 - C2 - C3 + C4)/4.0/dV/dS*2.0*Math.Sqrt(V0);
	
            // Volga
            dC2 = (V1 - 2.0*PricePDE + V2)/(dV*dV);
            double VolgaPDE = 4*Math.Sqrt(V0)*(dC2*Math.Sqrt(V0) + Vega1PDE/4.0/V0);

            // Theta
            T1 = IP.interp2(V,S,U,V0,S0);
            T2 = IP.interp2(V,S,u,V0,S0);
            double ThetaPDE = -(T1-T2)/dt;

            // Output the results
            Console.WriteLine("Stock price grid size  {0}", nS+1);
            Console.WriteLine("Volatility grid size   {0}",nV+1);
            Console.WriteLine("Number of time steps   {0}", nT);
            Console.WriteLine("---------------------------------------");
            Console.WriteLine("Greek             PDE        FiniteDiff");
            Console.WriteLine("---------------------------------------");
            Console.WriteLine("Price         {0,10:F5} {1,12:F5} ",PricePDE,PriceClosed);
            Console.WriteLine("Delta         {0,10:F5} {1,12:F5} ",DeltaPDE,DeltaFD);
            Console.WriteLine("Gamma         {0,10:F5} {1,12:F5} ",GammaPDE,GammaFD);
            Console.WriteLine("Vega #1       {0,10:F5} {1,12:F5} ",Vega1PDE,Vega1FD);
            Console.WriteLine("Vanna         {0,10:F5} {1,12:F5} ",VannaPDE,VannaFD);
            Console.WriteLine("Volga         {0,10:F5} {1,12:F5} ",VolgaPDE,VolgaFD);
            Console.WriteLine("Theta         {0,10:F5} {1,12:F5} ",ThetaPDE,ThetaFD);
            Console.WriteLine("---------------------------------------");
        }
    }
}



