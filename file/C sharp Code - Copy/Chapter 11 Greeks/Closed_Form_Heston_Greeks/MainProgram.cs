using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Closed_Form_Heston_Greeks
{
    class ClosedFormHestonGreeks
    {
        static void Main(string[] args)
        {
            // Classes
            HestonPrice HP = new HestonPrice();
            HestonGreeksAlgo HG = new HestonGreeksAlgo();
            BlackScholesPrice BS = new BlackScholesPrice();

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
            double kappa = 0.01;
            double theta = 0.07;
            double v0    = theta;
            double sigma = 1e-5;
            double rho = -0.8;
            double lambda = 0.0;

            double S = 100.0;
            double K = 100.0;
            double T = 0.15;
            double r = 0.05;
            double q = 0.03;
            string PutCall = "C";
            int trap = 1;

            OpSet settings = new OpSet();
            settings.S = S;
            settings.K = K;
            settings.T = T;
            settings.r = r;
            settings.q = q;
            settings.PutCall = PutCall;

            HParam param = new HParam();
            param.kappa = kappa;
            param.theta = theta;
            param.sigma = sigma;
            param.v0    = v0;
            param.rho   = rho;
            param.lambda = lambda;


            // Option price
            double HPrice  = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double vol  = Math.Sqrt(v0);
            double BSPrice = BS.BlackScholes(S,K,T,r,q,vol,PutCall);

            // Heston Greeks by closed form ======================================================================
            double Delta = HG.HestonGreeks(param,settings,X,W,"Delta");
            double Gamma = HG.HestonGreeks(param,settings,X,W,"Gamma");
            double Rho   = HG.HestonGreeks(param,settings,X,W,"Rho");
            double Theta = HG.HestonGreeks(param,settings,X,W,"Theta");
            double Vega1 = HG.HestonGreeks(param,settings,X,W,"Vega1");
            double Vega2 = HG.HestonGreeks(param,settings,X,W,"Vega2");
            double Vanna = HG.HestonGreeks(param,settings,X,W,"Vanna");
            double Volga = HG.HestonGreeks(param,settings,X,W,"Volga");

            // Heston Greeks by Finite Difference =================================================================
            // Delta
            double dS = 1.0;
            double D1 = HP.HestonPriceGaussLaguerre(S+dS,K,T,r,q,kappa,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double D2 = HP.HestonPriceGaussLaguerre(S-dS,K,T,r,q,kappa,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double DeltaFD = (D1-D2)/2.0/dS;

            // Gamma
            double GammaFD = (D1 - 2.0*HPrice + D2)/dS/dS;

            // Rho
            double dr = 1.0e-5;
            double R1 = HP.HestonPriceGaussLaguerre(S,K,T,r+dr,q,kappa,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double R2 = HP.HestonPriceGaussLaguerre(S,K,T,r-dr,q,kappa,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double RhoFD = (R1-R2)/2.0/dr;

            // Vega #1  
            double dv = 1e-5;
            double V1 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta,lambda,sigma,v0+dv,rho,trap,X,W,PutCall); ;
            double V2 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta,lambda,sigma,v0-dv,rho,trap,X,W,PutCall); ;
            double Vega1FD = (V1-V2)/2/dv*2.0*Math.Sqrt(v0);

            // Vega #2
            double dh = 1e-6;
            double H1 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta+dh,lambda,sigma,v0,rho,trap,X,W,PutCall); ;
            double H2 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta-dh,lambda,sigma,v0,rho,trap,X,W,PutCall); ;
            double Vega2FD = (H1-H2)/2.0/dh*2.0*Math.Sqrt(theta);

            // Theta
            double dt = 1.0e-5;
            double T1 = HP.HestonPriceGaussLaguerre(S,K,T+dt,r,q,kappa,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double T2 = HP.HestonPriceGaussLaguerre(S,K,T-dt,r,q,kappa,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double ThetaFD = -(T1-T2)/2.0/dt;

            // Vanna
            dv = 1e-5;
            dS = 1e-1;
            double C1 = HP.HestonPriceGaussLaguerre(S+dS,K,T,r,q,kappa,theta,lambda,sigma,v0+dv,rho,trap,X,W,PutCall);
            double C2 = HP.HestonPriceGaussLaguerre(S+dS,K,T,r,q,kappa,theta,lambda,sigma,v0-dv,rho,trap,X,W,PutCall);
            double C3 = HP.HestonPriceGaussLaguerre(S-dS,K,T,r,q,kappa,theta,lambda,sigma,v0+dv,rho,trap,X,W,PutCall);
            double C4 = HP.HestonPriceGaussLaguerre(S-dS,K,T,r,q,kappa,theta,lambda,sigma,v0-dv,rho,trap,X,W,PutCall);
            double VannaFD = (C1 - C2 - C3 + C4)/4.0/dv/dS*2.0*Math.Sqrt(v0);

            // Volga
            double dC2 = (V1 - 2*HPrice + V2)/(dv*dv);
            double VolgaFD = 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1FD/4.0/v0);

            // Corr
            double dp = 1e-4;
            double Co1 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta,lambda,sigma,v0,rho+dp,trap,X,W,PutCall);
            double Co2 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta,lambda,sigma,v0,rho-dp,trap,X,W,PutCall);
            double CorrFD = (Co1-Co2)/2.0/dp;

            // Volatility of variance
            double ds = 1e-2;
            double s1 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta,lambda,sigma+ds,v0,rho,trap,X,W,PutCall);
            double s2 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa,theta,lambda,sigma-ds,v0,rho,trap,X,W,PutCall);
            double SigmaFD = (s1-s2)/2.0/ds;

            // Mean reversion speed
            double dk = 1e-2;
            double K1 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa+dk,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double K2 = HP.HestonPriceGaussLaguerre(S,K,T,r,q,kappa-dk,theta,lambda,sigma,v0,rho,trap,X,W,PutCall);
            double KappaFD = (K1-K2)/2.0/dk;

            // Black Scholes Greeks ===================================================================================
            double BSDelta = BS.BSGreeks(PutCall,S,K,r,q,T,vol,"Delta");
            double BSGamma = BS.BSGreeks(PutCall,S,K,r,q,T,vol,"Gamma");
            double BSTheta = BS.BSGreeks(PutCall,S,K,r,q,T,vol,"Theta");
            double BSRho   = BS.BSGreeks(PutCall,S,K,r,q,T,vol,"Rho");
            double BSVega  = BS.BSGreeks(PutCall,S,K,r,q,T,vol,"Vega");
            double BSVanna = BS.BSGreeks(PutCall,S,K,r,q,T,vol,"Vanna");
            double BSVolga = BS.BSGreeks(PutCall,S,K,r,q,T,vol,"Volga");

            // Output the results
            Console.WriteLine("Exact Price             {0:F5}",HPrice);
            Console.WriteLine("Black-Scholes Price     {0:F5}",BSPrice);
            Console.WriteLine("-----------------------------------------------------");
            Console.WriteLine("Greek          Exact       Fin-Diff   Black-Scholes");
            Console.WriteLine("-----------------------------------------------------");
            Console.WriteLine("Delta    {0,12:F5} {1,12:F5} {2,12:F5}",Delta,DeltaFD,BSDelta);
            Console.WriteLine("Gamma    {0,12:F5} {1,12:F5} {2,12:F5}",Gamma,GammaFD,BSGamma);
            Console.WriteLine("Rho      {0,12:F5} {1,12:F5} {2,12:F5}",Rho,RhoFD,BSRho);
            Console.WriteLine("Theta    {0,12:F5} {1,12:F5} {2,12:F5}",Theta,ThetaFD,BSTheta);
            Console.WriteLine("Vega#1   {0,12:F5} {1,12:F5} {2,12:F5}",Vega1,Vega1FD,BSVega);
            Console.WriteLine("Vanna    {0,12:F5} {1,12:F5} {2,12:F5}",Vanna,VannaFD,BSVanna);
            Console.WriteLine("Volga    {0,12:F5} {1,12:F5} {2,12:F5}",Volga,VolgaFD,BSVolga);
            Console.WriteLine(" ");
            Console.WriteLine("Other Sensitivities");
            Console.WriteLine("-----------------------------------");
            Console.WriteLine("Vega#2   {0,12:F5} {1,12:F5}",Vega2,Vega2FD);
            Console.WriteLine("Corr     {0,25:E2}",CorrFD);
            Console.WriteLine("Sigma    {0,25:E2}",SigmaFD);
            Console.WriteLine("Kappa    {0,25:E2}",KappaFD);
            Console.WriteLine("-----------------------------------");
        }
    }
}

