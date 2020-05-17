using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Double_Heston_Greeks
{
    class DoubleHestonGreeks
    {
        static void Main(string[] args)
        {            // 32-point Gauss-Laguerre Abscissas and weights
            double[] X = new Double[32];
            double[] W = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    X[k] = double.Parse(bits[0]);
                    W[k] = double.Parse(bits[1]);
                }
            }
            DHGreeks DH = new DHGreeks();

            // Double Heston parameters
            DHParam param = new DHParam();
            double v01 = 0.6*0.6;
            double v02 = 0.7*0.7;
            param.v01 = v01;
            param.v02 = v02;
            param.sigma1 = 0.10;
            param.sigma2 = 0.20;
            param.kappa1 = 0.90;
            param.kappa2 = 1.20;
            param.rho1 = -0.5;
            param.rho2 = -0.5;
            param.theta1 = 0.10;
            param.theta2 = 0.15;

            // Option settings
            double S = 61.90;
            double r = 0.03;
            double T = 1;
            OpSet settings = new OpSet();
            settings.S = S;
            settings.K = S;
            settings.T = T;
            settings.r = r;
            settings.q = 0.0;
            settings.PutCall = "C";
            settings.trap = 1;

            // Calculate the price
            double Price = DH.DoubleHestonGreeks(param,settings,X,W,"Price");

            // Greeks by finite difference approximations
            // Delta and Gamma
            double ds = 0.1;
            double C0 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            settings.S = S + ds;
            double C1 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            settings.S = S - ds;
            double C1_ = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            double DeltaFD = (C1 - C1_)/2.0/ds;
            double GammaFD = (C1 - 2.0*C0 + C1_)/ds/ds;
            settings.S = S;

            // Rho
            double dr = 0.001;
            settings.r = r + dr;
            double R1  = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            settings.r = r - dr;
            double R1_ = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            double RhoFD = (R1 - R1_)/2.0/dr;
            settings.r = r;

            // Vega and Volga
            double dv = 0.001;
            param.v01 = v01 + dv;
            double V1  = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v01 = v01 - dv;
            double V1_ = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v01 = v01;
            double V0 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            double Vega11FD = (V1 - V1_)/2.0/dv*2.0*Math.Sqrt(v01);
            double dC1 = (V1 - 2*V0 + V1_)/dv/dv;
            double Volga1FD = 4.0*Math.Sqrt(v01)*(dC1*Math.Sqrt(v01) + Vega11FD/4.0/v01);

            param.v02 = v02 + dv;
            V1  = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v02 = v02 - dv;
            V1_ = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v02 = v02;
            V0 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            double Vega12FD = (V1 - V1_)/2.0/dv*2.0*Math.Sqrt(v02);
            dC1 = (V1 - 2*V0 + V1_)/dv/dv;
            double Volga2FD = 4.0*Math.Sqrt(v02)*(dC1*Math.Sqrt(v02) + Vega12FD/4.0/v02);

            // Theta
            double dt = 0.001;
            settings.T = T + dt;
            double T1  = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            settings.T = T - dt;
            double T1_ = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            double ThetaFD = -(T1 - T1_)/2.0/dt;
            settings.T = T;

            // Vanna
            param.v01 = v01+dv;
            settings.S = S+ds;
            double Va1 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v01 = v01+dv;
            settings.S = S-ds;
            double Va2 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v01 = v01-dv;
            settings.S = S+ds;
            double Va3 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v01 = v01-dv;
            settings.S = S-ds;
            double Va4 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v01 = v01;
            settings.S = S;
            double Vanna1FD = (Va1 - Va2 - Va3 + Va4)/4.0/dv/ds*2*Math.Sqrt(v01);

            param.v02 = v02+dv;
            settings.S = S+ds;
            Va1 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v02 = v02+dv;
            settings.S = S-ds;
            Va2 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v02 = v02-dv;
            settings.S = S+ds;
            Va3 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v02 = v02-dv;
            settings.S = S-ds;
            Va4 = DH.DoubleHestonGreeks(param,settings,X,W,"Price");
            param.v02 = v02;
            settings.S = S;
            double Vanna2FD = (Va1 - Va2 - Va3 + Va4)/4.0/dv/ds*2*Math.Sqrt(v02);

            // Greeks in analytic form
            double Delta = DH.DoubleHestonGreeks(param,settings,X,W,"Delta");
            double Gamma = DH.DoubleHestonGreeks(param,settings,X,W,"Gamma");
            double Rho   = DH.DoubleHestonGreeks(param,settings,X,W,"Rho");
            double Vega11 = DH.DoubleHestonGreeks(param,settings,X,W,"Vega11");
            double Vega12 = DH.DoubleHestonGreeks(param,settings,X,W,"Vega12");
            double Volga1 = DH.DoubleHestonGreeks(param,settings,X,W,"Volga1");
            double Volga2 = DH.DoubleHestonGreeks(param,settings,X,W,"Volga2");
            double Theta = DH.DoubleHestonGreeks(param,settings,X,W,"Theta");
            double Vanna1 = DH.DoubleHestonGreeks(param,settings,X,W,"Vanna1");
            double Vanna2 = DH.DoubleHestonGreeks(param,settings,X,W,"Vanna2");

            // Show that the double Heston PDE holds
            double Vega11p = Vega11/2.0/Math.Sqrt(v01);
            double Vega12p = Vega12/2.0/Math.Sqrt(v02);
            double Vanna1p = Vanna1/2.0/Math.Sqrt(v01);
            double Vanna2p = Vanna2/2.0/Math.Sqrt(v02);
            double Volga1p = (0.25/Math.Sqrt(v01)*Volga1 - Vega11p/2.0/Math.Sqrt(v01))/Math.Sqrt(v01);
            double Volga2p = (0.25/Math.Sqrt(v02)*Volga2 - Vega12p/2.0/Math.Sqrt(v02))/Math.Sqrt(v02);

            // The generator
            double A = (r-settings.q)*S*Delta 
                + param.kappa1*(param.theta1-v01)*Vega11p
                + param.kappa2*(param.theta2-v02)*Vega12p
                + 0.5*(v01+v02)*S*S*Gamma
                + param.rho1*param.sigma1*v01*S*Vanna1p
                + param.rho2*param.sigma2*v02*S*Vanna2p
                + 0.5*param.sigma1*param.sigma1*v01*Volga1p
                + 0.5*param.sigma2*param.sigma2*v02*Volga2p;
            double PDE = Theta + A - r*Price;
            
            // Output the results
            Console.WriteLine("------------------------------------------------");
            Console.WriteLine("Greek          Analytic   FiniteDiff    Error");
            Console.WriteLine("------------------------------------------------");
            Console.WriteLine("Price  {0,15:F4}",Price);
            Console.WriteLine("Delta  {0,15:F4} {1,10:F4} {2,12:E2}",Delta,DeltaFD,Math.Abs(Delta-DeltaFD));
            Console.WriteLine("Gamma  {0,15:F4} {1,10:F4} {2,12:E2}",Gamma,GammaFD,Math.Abs(Gamma-GammaFD));
            Console.WriteLine("Rho    {0,15:F4} {1,10:F4} {2,12:E2}",Rho,RhoFD,Math.Abs(Rho-RhoFD));
            Console.WriteLine("Theta  {0,15:F4} {1,10:F4} {2,12:E2}",Theta,ThetaFD,Math.Abs(Theta-ThetaFD));
            Console.WriteLine("Vega11 {0,15:F4} {1,10:F4} {2,12:E2}",Vega11,Vega11FD,Math.Abs(Vega11-Vega11FD));
            Console.WriteLine("Vega12 {0,15:F4} {1,10:F4} {2,12:E2}",Vega12,Vega12FD,Math.Abs(Vega12-Vega12FD));
            Console.WriteLine("Vanna1 {0,15:F4} {1,10:F4} {2,12:E2}",Vanna1,Vanna1FD,Math.Abs(Vanna1-Vanna1FD));
            Console.WriteLine("Vanna2 {0,15:F4} {1,10:F4} {2,12:E2}",Vanna2,Vanna2FD,Math.Abs(Vanna2-Vanna2FD));
            Console.WriteLine("Volga1 {0,15:F4} {1,10:F4} {2,12:E2}",Volga1,Volga1FD,Math.Abs(Volga1-Volga1FD));
            Console.WriteLine("Volga2 {0,15:F4} {1,10:F4} {2,12:E2}",Volga2,Volga2FD,Math.Abs(Volga2-Volga2FD));
            Console.WriteLine("------------------------------------------------");
            Console.WriteLine("PDE Value   {0,10:E2}",PDE);
            Console.WriteLine("------------------------------------------------");
        }
    }
}
