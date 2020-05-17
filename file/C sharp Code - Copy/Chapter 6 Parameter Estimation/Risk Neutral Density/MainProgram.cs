using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Risk_Neutral_Density
{
    class RiskNeutralDensity
    {
        static void Main(string[] args)
        {
            // Classes
            RND RND = new RND();
            HestonPriceMD HPMD = new HestonPriceMD();

            // 32-point Gauss-Laguerre Abscissas and weights
            double[] XGLe = new Double[32];
            double[] WGLe = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLegendre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    XGLe[k] = double.Parse(bits[0]);
                    WGLe[k] = double.Parse(bits[1]);
                }
            }

            // Maturities and strikes
            double[] T = new double[] { 45.0/365.0,98.0/365.0,261.0/365.0,348.0/365.0 };
            int NT = T.Length;

            // Option price settings
            OpSet settings = new OpSet();
            settings.S = 137.14;
            settings.r = 0.0010;
            settings.q = 0.0068;
            settings.trap = 1;
            settings.PutCall = "C";

            // Heston RMSE parameters from SP data, Table 6.1
            HParam param = new HParam();
            param.kappa =  8.9957;
            param.theta =  0.0567;
            param.sigma =  2.0000;
            param.v0    =  0.0415;
            param.rho   = -0.7804;
            param.lambda = 0;

            // Potential integration domain 
            double lo = 1e-10;
            double hi = 500;
            int Ndomains = 50;
            double tol = 1e-6;
            double dA = (hi-lo)/Convert.ToDouble(Ndomains);
            double[] A = new double[Ndomains];
            for(int j=0;j<=Ndomains-1;j++)
                A[j] = lo + j*dA;

            // Initialize the domains and number of points
            double[] Domain = new double[NT];
            double[] Points = new double[NT];
            double[] Area = new double[NT];

            // First maturity =====================================================================================================
            int NK1 = 100;
            settings.T = T[0];
            double LoK = 90.0;
            double HiK = 170.0;
            double dK = (HiK - LoK)/NK1;
            double[] K1 = new double[NK1+1];
            double[] CallPrice1 = new double[NK1+1];
            double[] points = new double[NK1+1];
            double[] Upper = new double[NK1+1];
            OutputMD MDoutput;
            for(int k=0;k<=NK1;k++)
            {
                K1[k] = LoK + k*dK;
                settings.K = K1[k];
                MDoutput = HPMD.HestonPriceGaussLegendreMD(param,settings,XGLe,WGLe,A,tol);
                CallPrice1[k] = MDoutput.Price;
                points[k] = MDoutput.Npoints;
                Upper[k] = MDoutput.upper;
            }

            fS RNDoutput = RND.ExtractRND(K1,CallPrice1);
            double[] RND1 = RNDoutput.RND;
            double[] Strike1 = RNDoutput.K;
            double dh = Strike1[1] - Strike1[0];
            Domain[0] = Upper.Max();
            Points[0] = points.Max();
            Area[0] = RND1.Sum() * dh;

            // Second maturity =====================================================================================================
            int NK2 = 100;
            settings.T = T[1];
            LoK = 80.0;
            HiK = 180.0;
            dK = (HiK - LoK)/NK2;
            double[] K2 = new double[NK2+1];
            double[] CallPrice2 = new double[NK2+1];
            points = new double[NK2+1];
            Upper = new double[NK2+1];
            for(int k=0;k<=NK2;k++)
            {
                K2[k] = LoK + k*dK;
                settings.K = K2[k];
                MDoutput = HPMD.HestonPriceGaussLegendreMD(param,settings,XGLe,WGLe,A,tol);
                CallPrice2[k] = MDoutput.Price;
                points[k] = MDoutput.Npoints;
                Upper[k] = MDoutput.upper;
            }
            RNDoutput = RND.ExtractRND(K2,CallPrice2);
            double[] RND2 = RNDoutput.RND;
            double[] Strike2 = RNDoutput.K;
            dh = Strike2[1] - Strike2[0];
            Domain[1] = Upper.Max();
            Points[1] = points.Max();
            Area[1] = RND2.Sum() * dh;

            // Third maturity =====================================================================================================
            int NK3 = 150;
            settings.T = T[2];
            LoK = 40.0;
            HiK = 220.0;
            dK = (HiK - LoK)/NK3;
            double[] K3 = new double[NK3+1];
            double[] CallPrice3 = new double[NK3+1];
            points = new double[NK3+1];
            Upper = new double[NK3+1];
            for(int k=0;k<=NK3;k++)
            {
                K3[k] = LoK + k*dK;
                settings.K = K3[k];
                MDoutput = HPMD.HestonPriceGaussLegendreMD(param,settings,XGLe,WGLe,A,tol);
                CallPrice3[k] = MDoutput.Price;
                points[k] = MDoutput.Npoints;
                Upper[k] = MDoutput.upper;
            }
            RNDoutput = RND.ExtractRND(K3,CallPrice3);
            double[] RND3 = RNDoutput.RND;
            double[] Strike3 = RNDoutput.K;
            dh = Strike3[1] - Strike3[0];
            Domain[2] = Upper.Max();
            Points[2] = points.Max();
            Area[2] = RND3.Sum() * dh;

            // Fourth maturity =====================================================================================================
            int NK4 = 200;
            settings.T = T[3];
            LoK = 20.0;
            HiK = 240.0;
            dK = (HiK - LoK)/NK4;
            double[] K4 = new double[NK4+1];
            double[] CallPrice4 = new double[NK4+1];
            points = new double[NK4+1];
            Upper = new double[NK4+1];
            for(int k=0;k<=NK4;k++)
            {
                K4[k] = LoK + k*dK;
                settings.K = K4[k];
                MDoutput = HPMD.HestonPriceGaussLegendreMD(param,settings,XGLe,WGLe,A,tol);
                CallPrice4[k] = MDoutput.Price;
                points[k] = MDoutput.Npoints;
                Upper[k] = MDoutput.upper;
            }
            RNDoutput = RND.ExtractRND(K4,CallPrice4);
            double[] RND4 = RNDoutput.RND;
            double[] Strike4 = RNDoutput.K;
            dh = Strike4[1] - Strike4[0];
            Domain[3] = Upper.Max();
            Points[3] = points.Max();
            Area[3] = RND4.Sum() * dh;

            // Output the results
            Console.WriteLine("Maturity    Area      UpperLimit  Points");
            Console.WriteLine("----------------------------------------");
            for(int t=0;t<=NT-1;t++)
                Console.WriteLine(" {0,3:F0} {1,12:F5} {2,10:F0} {3,10:F0}",T[t]*365.0,Area[t],Domain[t],Points[t]);
            Console.WriteLine("----------------------------------------");
        }
    }
}
