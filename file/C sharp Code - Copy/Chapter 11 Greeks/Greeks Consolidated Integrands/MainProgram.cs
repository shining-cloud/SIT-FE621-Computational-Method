using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Greeks_Consolidated_Integrands
{
    static class ConsolidatedGreeks
    {
        static void Main(string[] args)
        {
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
            // Option settings
            double S = 100.0;
            double K = 103.0;
            double r = 0.05;
            double q = 0.03;
            double T = 0.5;
            int trap = 1;

            // Heston parameters
            double rho = -0.8;
            double kappa  = 5;
            double lambda = 0;
            double sigma  = 0.5;
            double v0 = 0.07;
            double theta = 0.05;

            // Greeks by finite diferences
            double Price = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Price");

            double ds = 0.1;
            double C1  = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S+ds,r,q,v0,trap,x,w,"Price");
            double C0  = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Price");
            double C1_ = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S-ds,r,q,v0,trap,x,w,"Price");
            double DeltaFD = (C1 - C1_)/2.0/ds;
            double GammaFD = (C1 - 2.0*C0 + C1_)/ds/ds;

            double dt = 0.01;
            double T1  = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T+dt,K,S,r,q,v0,trap,x,w,"Price");
            double T1_ = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T-dt,K,S,r,q,v0,trap,x,w,"Price");
            double ThetaFD = -(T1 - T1_)/2.0/dt;

            double dr = 0.01;
            double R1  = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r+dr,q,v0,trap,x,w,"Price");
            double R1_ = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r-dr,q,v0,trap,x,w,"Price");
            double RhoFD = (R1 - R1_)/2.0/dr;

            double dv = 1e-3;
            double V1  = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0+dv,trap,x,w,"Price");
            double V0  = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Price");
            double V1_ = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0-dv,trap,x,w,"Price");
            double Vega1FD = (V1 - V1_)/2.0/dv*2.0*Math.Sqrt(v0);
            double dC2 = (V1 - 2.0*V0 + V1_)/dv/dv;
            double VolgaFD = 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1FD/4.0/v0);

            double Va1 = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S+ds,r,q,v0+dv,trap,x,w,"Price");
            double Va2 = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S+ds,r,q,v0-dv,trap,x,w,"Price");
            double Va3 = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S-ds,r,q,v0+dv,trap,x,w,"Price");
            double Va4 = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S-ds,r,q,v0-dv,trap,x,w,"Price");
            double VannaFD = (Va1 - Va2 - Va3 + Va4)/4.0/dv/ds*2.0*Math.Sqrt(v0);

            // Greeks in Analytic form
            double Delta = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Delta");
            double Gamma = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Gamma");
            double Theta = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Theta");
            double Rho   = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Rho");
            double Vega1 = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Vega1");
            double Vanna = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Vanna");
            double Volga = HG.HestonGreekConsolidated(kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,x,w,"Volga");

            // Output the results
            Console.WriteLine("Greek           Analytic      Finite Diff");
            Console.WriteLine("-----------------------------------------");
            Console.WriteLine("Price   {0,15:F5}",Price);
            Console.WriteLine("Delta   {0,15:F5} {1,15:F5}",Delta,DeltaFD);
            Console.WriteLine("Gamma   {0,15:F5} {1,15:F5}",Gamma,GammaFD);
            Console.WriteLine("Rho     {0,15:F5} {1,15:F5}",Rho,RhoFD);
            Console.WriteLine("Theta   {0,15:F5} {1,15:F5}",Theta,ThetaFD);
            Console.WriteLine("Vega1   {0,15:F5} {1,15:F5}",Vega1,Vega1FD);
            Console.WriteLine("Vanna   {0,15:F5} {1,15:F5}",Vanna,VannaFD);
            Console.WriteLine("Volga   {0,15:F5} {1,15:F5}",Volga,VolgaFD);
            Console.WriteLine("-----------------------------------------");



        }
    }
}


