using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Carr_and_Madan_Greeks
{
    class CarrMadanGreeks
    {
        static void Main(string[] args)
        {
            HestonPrice HP = new HestonPrice();
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
            double S = 100.0;				    // Spot Price
            double K = 105.0;				    // Strike Price
            double T = 0.5;			            // Maturity in Years
            double r = 0.05;					// Interest Rate
            double q = 0.0;                     // Dividend yield
            double rho = -0.7;					// Heston Parameter: Correlation
            double kappa = 2.0;					// Heston Parameter 
            double theta = 0.06;				// Heston Parameter 
            double lambda = 0.0;				// Heston Parameter 
            double sigma = 0.1;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.06;					// Heston Parameter: Current Variance
            int trap = 1;                       // 1="Little Trap" characteristic function
            double alpha = 1.75;                // Carr Madan damping factor
            string PutCall = "C";

            // Calculate the Greeks by finite differences
            double PriceFD = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);

            // Delta
            double ds = 0.1;
            double S1  = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S+ds,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double S0  = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double S1_ = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S-ds,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double DeltaFD = (S1 - S1_)/2.0/ds;
            double GammaFD = (S1 - 2.0*S0 + S1_)/ds/ds;

            // Theta
            double dt = 1e-4;
            double T1  = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r,q,T+dt,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double T1_ = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r,q,T-dt,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double ThetaFD = (T1 - T1_)/2.0/dt;

            // Rho
            double dr = 1e-5;
            double R1  = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r+dr,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double R1_ = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r-dr,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double RhoFD = (R1 - R1_)/2.0/dr;

            // Vega1
            double dv0 = 1e-4;
            double Ve1  = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r,q,T,kappa,theta,sigma,v0+dv0,lambda,rho,x,w,trap);
            double Ve1_ = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r,q,T,kappa,theta,sigma,v0-dv0,lambda,rho,x,w,trap);
            double Vega1FD = (Ve1 - Ve1_)/2.0/dv0*2.0*Math.Sqrt(v0);

            // Vanna
            double Va1 = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S+ds,K,r,q,T,kappa,theta,sigma,v0+dv0,lambda,rho,x,w,trap);
            double Va2 = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S+ds,K,r,q,T,kappa,theta,sigma,v0-dv0,lambda,rho,x,w,trap);
            double Va3 = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S-ds,K,r,q,T,kappa,theta,sigma,v0+dv0,lambda,rho,x,w,trap);
            double Va4 = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S-ds,K,r,q,T,kappa,theta,sigma,v0-dv0,lambda,rho,x,w,trap);
            double VannaFD = (Va1 - Va2 - Va3 + Va4)/4/dv0/ds*2.0*Math.Sqrt(v0);

            // Volga
            double Ve0 = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,alpha,S,K,r,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double dC2 = (Ve1 - 2.0*Ve0 + Ve1_)/dv0/dv0;
            double VolgaFD = 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1FD/4.0/v0);

            // Greeks by closed form
            double Price = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Price",x,w);
            double Delta = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Delta",x,w);
            double Gamma = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Gamma",x,w);
            double Theta = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Theta",x,w);
            double Rho   = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Rho",x,w);
            double Vega1 = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Vega1",x,w);
            double Vanna = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Vanna",x,w);
            double Volga = HG.CarrMadanGreeks(alpha,kappa,theta,lambda,rho,sigma,T,K,S,r,v0,trap,"Volga",x,w);

            // Output the results
            Console.WriteLine("Carr-Madan Greeks by finite differences and by closed form");
            Console.WriteLine("-----------------------------------");
            Console.WriteLine("Greek          F.D.        Closed");
            Console.WriteLine("-----------------------------------");
            Console.WriteLine("Call     {0,10:F4}    {1,10:F4}",PriceFD,Price);
            Console.WriteLine("Delta    {0,10:F4}    {1,10:F4}",DeltaFD,Delta);
            Console.WriteLine("Gamma    {0,10:F4}    {1,10:F4}",GammaFD,Gamma);
            Console.WriteLine("Theta    {0,10:F4}    {1,10:F4}",ThetaFD,Theta);
            Console.WriteLine("Rho      {0,10:F4}    {1,10:F4}",RhoFD,Rho);
            Console.WriteLine("Vega1    {0,10:F4}    {1,10:F4}",Vega1FD,Vega1);
            Console.WriteLine("Vanna    {0,10:F4}    {1,10:F4}",VannaFD,Vanna);
            Console.WriteLine("Volga    {0,10:F4}    {1,10:F4}",VolgaFD,Volga);
            Console.WriteLine("-----------------------------------");
        }
    }
}
