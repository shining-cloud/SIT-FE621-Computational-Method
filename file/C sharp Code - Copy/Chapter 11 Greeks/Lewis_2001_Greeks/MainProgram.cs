using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_2001_Greeks
{
    class Lewis2001Greeks
    {
        static void Main(string[] args)
        {
            LewisGreeks LG = new LewisGreeks();
            
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
            double S = 100.0;				    // Spot Price
            double K = 103.0;                   // Strike Price
            double T = 0.5; 			        // Maturity in Years
            double r = 0.05;					// Interest Rate
            double q = 0.03;                    // Dividend yield
            double rho = -0.8;					// Heston Parameter: Correlation
            double kappa = 5;					// Heston Parameter 
            double theta = 0.05;				// Heston Parameter 
            double sigma = 0.5;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.07;					// Heston Parameter: Current Variance
            int trap = 1;                       // 1="Little Trap" characteristic function

            // Greeks by finite difference
            double Price = LG.LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,"Price");

            // Delta and Gamma
            double ds = 0.1;
            double C1  = LG.LewisGreeks311(S+ds,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double C0  = LG.LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double C1_ = LG.LewisGreeks311(S-ds,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double DeltaFD = (C1 - C1_)/2.0/ds;
            double GammaFD = (C1 - 2.0*C0 + C1_)/ds/ds;

            // Rho
            double dr = 1e-4;
            double R1  = LG.LewisGreeks311(S,K,r+dr,q,T,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double R1_ = LG.LewisGreeks311(S,K,r-dr,q,T,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double RhoFD = (R1 - R1_)/2.0/dr;

            // Theta
            double dt = 1e-4;
            double T1  = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double T1_ = LG.LewisGreeks311(S,K,r,q,T-dt,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double ThetaFD = -(T1 - T1_)/2.0/dt;

            // Vega and Volga
            double dv = 1e-4;
            double V1  = LG.LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0+dv,trap,x,w,"Price");
            double V0  = LG.LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,"Price");
            double V1_ = LG.LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0-dv,trap,x,w,"Price");
            double Vega1FD = (V1 - V1_)/2.0/dv*2.0*Math.Sqrt(v0);
            double dC2 = (V1 - 2*V0 + V1_)/dv/dv;
            double VolgaFD = 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1FD/4.0/v0);

            // Vanna
            double Va1 = LG.LewisGreeks311(S+ds,K,r,q,T,theta,kappa,sigma,rho,v0+dv,trap,x,w,"Price");
            double Va2 = LG.LewisGreeks311(S+ds,K,r,q,T,theta,kappa,sigma,rho,v0-dv,trap,x,w,"Price");
            double Va3 = LG.LewisGreeks311(S-ds,K,r,q,T,theta,kappa,sigma,rho,v0+dv,trap,x,w,"Price");
            double Va4 = LG.LewisGreeks311(S-ds,K,r,q,T,theta,kappa,sigma,rho,v0-dv,trap,x,w,"Price");
            double VannaFD = (Va1 - Va2 - Va3 + Va4)/4.0/dv/ds*2.0*Math.Sqrt(v0);

            // Closed form Greeks
            double Delta = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Delta");
            double Gamma = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Gamma");
            double Rho = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Rho");
            double Theta = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Theta");
            double Vega1 = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Vega1");
            double Volga = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Volga");
            double Vanna = LG.LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,"Vanna");

            Console.WriteLine("Greek       Analytic   Finite Diff");
            Console.WriteLine("----------------------------------");
            Console.WriteLine("Price   {0,10:F4}",Price);
            Console.WriteLine("Delta   {0,10:F4} {1,12:F4}",Delta,DeltaFD);
            Console.WriteLine("Gamma   {0,10:F4} {1,12:F4}",Gamma,GammaFD);
            Console.WriteLine("Rho     {0,10:F4} {1,12:F4}",Rho,RhoFD);
            Console.WriteLine("Theta   {0,10:F4} {1,12:F4}",Theta,ThetaFD);
            Console.WriteLine("Vega1   {0,10:F4} {1,12:F4}",Vega1,Vega1FD);
            Console.WriteLine("Vanna   {0,10:F4} {1,12:F4}",Vanna,VannaFD);
            Console.WriteLine("Volga   {0,10:F4} {1,12:F4}",Volga,VolgaFD);
            Console.WriteLine("----------------------------------");


        }
    }
}
