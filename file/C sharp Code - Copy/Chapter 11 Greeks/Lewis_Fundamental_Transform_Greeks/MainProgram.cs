using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Lewis_Fundamental_Transform_Greeks
{
    class Lewis2000Greeks
    {
        static void Main(string[] args)
        {
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
            double K = 100.0;				    // Strike Price
            double T = 0.5;			            // Maturity in Years
            double r = 0.05;					// Interest Rate
            double q = 0.03;                    // Dividend yield
            double rho = -0.8;					// Heston Parameter: Correlation
            double kappa = 5;					// Heston Parameter 
            double theta = 0.05;				// Heston Parameter 
            double sigma = 0.5;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.07;					// Heston Parameter: Current Variance

            // First expression C1(K) for the Lewis Greeks =================================================
            double ki1 = 1.5;
            int form1 = 1;

            // Price
            LewisGreeks LG = new LewisGreeks();
            double PriceFD1 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");

            // Delta and gamma
            double ds = 1;
            double C1  = LG.HestonLewisGreekPrice(S+ds,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double C0  = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double C1_ = LG.HestonLewisGreekPrice(S-ds,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double DeltaFD1 = (C1 - C1_)/2.0/ds;
            double GammaFD1 = (C1 - 2.0*C0 + C1_)/ds/ds;

            // Rho
            double dr = 1e-5;
            double R1  = LG.HestonLewisGreekPrice(S,K,r+dr,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double R1_ = LG.HestonLewisGreekPrice(S,K,r-dr,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double RhoFD1 = (R1 - R1_)/2.0/dr;

            // Theta
            double dt = 1e-5;
            double T1  = LG.HestonLewisGreekPrice(S,K,r,q,v0,T+dt,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double T1_ = LG.HestonLewisGreekPrice(S,K,r,q,v0,T-dt,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double ThetaFD1 = -(T1 - T1_)/2.0/dt;

            // Vega and Volga
            double dv = 1e-4;
            double Ve1  = LG.HestonLewisGreekPrice(S,K,r,q,v0+dv,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double Ve0  = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double Ve1_ = LG.HestonLewisGreekPrice(S,K,r,q,v0-dv,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double Vega1FD1 = (Ve1 - Ve1_)/2.0/dv*2.0*Math.Sqrt(v0);
            double dC2 = (Ve1 - 2*Ve0 + Ve1_)/(dv*dv);
            double VolgaFD1 = 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1FD1/4.0/v0);

            // Vanna
            double Va1  = LG.HestonLewisGreekPrice(S+ds,K,r,q,v0+dv,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double Va2  = LG.HestonLewisGreekPrice(S+ds,K,r,q,v0-dv,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double Va3  = LG.HestonLewisGreekPrice(S-ds,K,r,q,v0+dv,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double Va4  = LG.HestonLewisGreekPrice(S-ds,K,r,q,v0-dv,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double VannaFD1 = (Va1 - Va2 - Va3 + Va4)/4.0/dv/ds*2.0*Math.Sqrt(v0);

            // Closed form Greeks for C1(K)
            double Price1 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Price");
            double Delta1 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Delta");
            double Gamma1 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Gamma");
            double Rho1   = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Rho");
            double Theta1 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Theta");
            double Vega11 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Vega1");
            double Vanna1 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Vanna");
            double Volga1 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki1,theta,kappa,sigma,rho,form1,x,w,"Volga");


            // Second expression C2(K) for the Lewis Greeks =================================================
            double ki2 = 0.5;
            int form2 = 2;

            // Price
            double PriceFD2 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");

            // Delta and gamma
            C1  = LG.HestonLewisGreekPrice(S+ds,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            C0  = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            C1_ = LG.HestonLewisGreekPrice(S-ds,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            double DeltaFD2 = (C1 - C1_)/2.0/ds;
            double GammaFD2 = (C1 - 2.0*C0 + C1_)/ds/ds;

            // Rho
            R1  = LG.HestonLewisGreekPrice(S,K,r+dr,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            R1_ = LG.HestonLewisGreekPrice(S,K,r-dr,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            double RhoFD2 = (R1 - R1_)/2.0/dr;

            // Theta
            T1  = LG.HestonLewisGreekPrice(S,K,r,q,v0,T+dt,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            T1_ = LG.HestonLewisGreekPrice(S,K,r,q,v0,T-dt,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            double ThetaFD2 = -(T1 - T1_)/2.0/dt;

            // Vega and Volga
            Ve1  = LG.HestonLewisGreekPrice(S,K,r,q,v0+dv,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            Ve0  = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            Ve1_ = LG.HestonLewisGreekPrice(S,K,r,q,v0-dv,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            double Vega1FD2 = (Ve1 - Ve1_)/2.0/dv*2.0*Math.Sqrt(v0);
            dC2 = (Ve1 - 2*Ve0 + Ve1_)/(dv*dv);
            double VolgaFD2 = 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1FD1/4.0/v0);

            // Vanna
            Va1  = LG.HestonLewisGreekPrice(S+ds,K,r,q,v0+dv,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            Va2  = LG.HestonLewisGreekPrice(S+ds,K,r,q,v0-dv,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            Va3  = LG.HestonLewisGreekPrice(S-ds,K,r,q,v0+dv,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            Va4  = LG.HestonLewisGreekPrice(S-ds,K,r,q,v0-dv,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            double VannaFD2 = (Va1 - Va2 - Va3 + Va4)/4.0/dv/ds*2.0*Math.Sqrt(v0);

            // Closed form Greeks for C2(K)
            double Price2 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Price");
            double Delta2 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Delta");
            double Gamma2 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Gamma");
            double Rho2   = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Rho");
            double Theta2 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Theta");
            double Vega12 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Vega1");
            double Vanna2 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Vanna");
            double Volga2 = LG.HestonLewisGreekPrice(S,K,r,q,v0,T,ki2,theta,kappa,sigma,rho,form2,x,w,"Volga");

            // Output the results
            Console.WriteLine("Comparison of Lewis (2000) C1(K) and C2(K) Greeks");
            Console.WriteLine(" ");
            Console.WriteLine("          ------- C1(K) -------      ------- C2(K) -------");
            Console.WriteLine("Greek     Analytic  Finite Diff      Analytic  Finite Diff");
            Console.WriteLine("-----------------------------------------------------------");
            Console.WriteLine("Price  {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Price1,PriceFD1,Price2,PriceFD2);
            Console.WriteLine("Delta  {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Delta1,DeltaFD1,Delta2,DeltaFD2);
            Console.WriteLine("Gamma  {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Gamma1,GammaFD1,Gamma2,GammaFD2);
            Console.WriteLine("Rho    {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Rho1,  RhoFD1  ,Rho2,  RhoFD2  );
            Console.WriteLine("Theta  {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Theta1,ThetaFD1,Theta2,ThetaFD2);
            Console.WriteLine("Vega1  {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Vega11,Vega1FD1,Vega12,Vega1FD2);
            Console.WriteLine("Vanna  {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Vanna1,VannaFD1,Vanna2,VannaFD2);
            Console.WriteLine("Volga  {0,10:F4} {1,10:F4} {2,15:F4} {3,10:F4}",Volga1,VolgaFD1,Volga2,VolgaFD2);
            Console.WriteLine("-----------------------------------------------------------");
        }
    }
}
