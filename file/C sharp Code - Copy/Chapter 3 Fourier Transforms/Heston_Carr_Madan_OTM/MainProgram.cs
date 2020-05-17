using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Carr_Madan_OTM
{
    class HestonCM_OTM
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
            double S = 1.0;				        // Spot Price
            double T = 1.0 ;			        // Maturity in Years
            double r = 0.03;					// Interest Rate
            double kappa = 2.0;					// Heston Parameter 
            double theta = 0.25;				// Heston Parameter 
            double sigma = 0.3;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.05;					// Heston Parameter: Current Variance
            double rho = -0.8;					// Heston Parameter: Correlation
            double lambda = 0.0;				// Heston Parameter 
            int trap = 1;                       // 1="Little Trap" characteristic function
            double alpha = 1.1;                 // Damping factor

            // Calculate the OTM Puts
            HestonPrice HP = new HestonPrice();
            string PutCall = "P";
            double Kp = 0.95;
            double HestonPut = HP.HestonPriceGaussLaguerre("Heston",PutCall,S,Kp,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);
            double CMPut     = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,S,Kp,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);
            double CMPutDamp = HP.HestonPriceGaussLaguerre("CarrMadanDamped",PutCall,S,Kp,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);

            // Calculate the OTM Calls
            PutCall = "C";
            double Kc = 1.05;
            double HestonCall = HP.HestonPriceGaussLaguerre("Heston",PutCall,S,Kc,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);
            double CMCall     = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,S,Kc,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);
            double CMCallDamp = HP.HestonPriceGaussLaguerre("CarrMadanDamped",PutCall,S,Kc,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);
            
            // Output the results
            Console.WriteLine("Spot = {0:F2}",S);
            Console.WriteLine("Method                OTM Put K = {0,4:F2}      OTM Call = {1,4:F2}",Kp,Kc);
            Console.WriteLine("-------------------------------------------------------------");
            Console.WriteLine("Heston                {0,10:F4} {1,20:F4}",HestonPut,HestonCall);
            Console.WriteLine("Carr Madan Undamped   {0,10:F4} {1,20:F4}",CMPut,CMCall);
            Console.WriteLine("Carr Madan Damped     {0,10:F4} {1,20:F4}",CMPutDamp,CMCallDamp);
            Console.WriteLine(" ");

            // Define the spot and strikes
            S = 25.0;

            // Calculate OTM puts
            PutCall = "P";
            Kp = 20.0;
            HestonPut = HP.HestonPriceGaussLaguerre("Heston",PutCall,S,Kp,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);
            CMPut     = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,S,Kp,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);

            // Calculate the OTM Calls
            PutCall = "C";
            Kc = 30.0;
            HestonCall = HP.HestonPriceGaussLaguerre("Heston",PutCall,S,Kc,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);
            CMCall     = HP.HestonPriceGaussLaguerre("CarrMadan",PutCall,S,Kc,r,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap,alpha);

            // Output the results
            Console.WriteLine("Spot = {0:F0}",S);
            Console.WriteLine("Method                OTM Put K = {0,4:F2}      OTM Call = {1,4:F2}",Kp,Kc);
            Console.WriteLine("-------------------------------------------------------------");
            Console.WriteLine("Heston                {0,10:F4} {1,20:F4}",HestonPut,HestonCall);
            Console.WriteLine("Carr Madan Undamped   {0,10:F4} {1,20:F4}",CMPut,CMCall);
            Console.WriteLine(" ");
        }
    }
}

