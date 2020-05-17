using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Beliaeva_Nawalkha_Tree
{
    partial class BNTree
    {
        static void Main(string[] args)
        {
            // Input the Gauss Laguerre abscissas and weights
            float[] X = new float[32];
            float[] W = new float[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
            {
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    X[k] = float.Parse(bits[0]);
                    W[k] = float.Parse(bits[1]);
                }
            }

            // Clarke and Parott (2010) example  ------------------------------------------------------------------------------------
            // Heston parameters
            float kappa = 5.0f;
            float theta = 0.16f;
            float sigma = 0.9f;
            float v0 = 0.0625f;
            float rho = 0.1f;

            // Heston parameters in a structure
            HParam param;
            param.kappa = kappa;
            param.theta = theta;
            param.sigma = sigma;
            param.v0 = v0;
            param.rho = rho;
            int trap = 1;

            // Option features
            string PutCall = "P";
            float Strike = 10.0f;

            // Clarke and Parrott (1999) settings with prices from Ikonen and Tiovanen (2008)
            float T = 0.25f;
            float rf = 0.1f;
            float q = 0.0f;
            float[] S0 = new float[5] { 8.0f,9.0f,10.0f,11.0f,12.0f };
            float[] TruePrice = new float[5] { 2.0f,1.107641f,0.520030f,0.213668f,0.082036f };

            // Tree settings
            int NT = 50;
            float threshold = 0.01f;

            // Arrays and structures for the American and European prices
            BivariateStruct output = new BivariateStruct();
            float[] AmerPrice = new float[5];
            float[] EuroPrice = new float[5];
            float[] CVAmerPrice = new float[5];
            float[] EuroClosedPrice = new float[5];

            // Create the American and European prices from the tree
            MainTree MT = new MainTree();
            for(int s=0;s<=4;s++)
            {
                output = MT.BuildBivariateTree(S0[s],PutCall,Strike,T,rf,NT,kappa,theta,sigma,v0,rho,threshold);
                AmerPrice[s] = output.AmerPrice;
                EuroPrice[s] = output.EuroPrice;

                // Closed Form European price
                EuroClosedPrice[s] = HestonPriceGaussLaguerre(param,S0[s],Strike,rf,q,T,trap,PutCall,X,W);

                // American price by control variate method
                CVAmerPrice[s] = AmerPrice[s] + (EuroClosedPrice[s] - EuroPrice[s]);
            }

            // Output the results
            Console.WriteLine("---------------------------------------------------------------------------");
            Console.WriteLine("Heston American Option using Believa-Nawalkha tree -- Clarke Parrott ");
            Console.WriteLine("Using {0:0} time steps",NT);
            Console.WriteLine("---------------------------------------------------------------------------");
            Console.WriteLine(" Spot  AmerPutTree  AmerPutCV  TruePrice");
            Console.WriteLine("--------------------------------------------");
            for(int s=0;s<=4;s++)
                Console.WriteLine(" {0,3:0} {1,10:F6}  {2,10:F6} {3,10:F6}",S0[s],AmerPrice[s],CVAmerPrice[s],TruePrice[s]);
            Console.WriteLine("--------------------------------------------");
        }
    }
}


/*
            // Belieava and Nawalkha (2010) Exhibit 14 ------------------------------------------------------------------------------------
 
            // Heston parameters
            float kappa = 3.0f;
            float theta = 0.04f;
            float sigma = 0.1f;

            // Heston parameters in a structure
            HParam param;
            param.kappa = kappa;
            param.theta = theta;
            param.sigma = sigma;
            int trap = 1;

            // Option features
            float rf = 0.05f;
            float q = 0.0f;
            string PutCall = "P";
            float Strike = 100.0f;

            // Beliaeva and Nawalkha (2010) settings
            float[] S0 = new float[36] { 
                90.0f,100.0f,110.0f,90.0f,100.0f,110.0f,90.0f,100.0f,110.0f,
                90.0f,100.0f,110.0f,90.0f,100.0f,110.0f,90.0f,100.0f,110.0f,
                90.0f,100.0f,110.0f,90.0f,100.0f,110.0f,90.0f,100.0f,110.0f,
                90.0f,100.0f,110.0f,90.0f,100.0f,110.0f,90.0f,100.0f,110.0f };
            float[] T = new float[36] { 
                0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,0.0833f,
                0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,0.2500f,
                0.5000f,0.5000f,0.5000f,0.5000f,0.5000f,0.5000f,0.5000f,0.5000f,0.5000f,0.5000f,0.5000f,0.5000f };
            float[] rho = new float[36] { 
                -0.1000f,-0.1000f,-0.1000f,-0.7000f,-0.7000f,-0.7000f,-0.1000f,-0.1000f,-0.1000f,-0.7000f,-0.7000f,-0.7000f,
                -0.1000f,-0.1000f,-0.1000f,-0.7000f,-0.7000f,-0.7000f,-0.1000f,-0.1000f,-0.1000f,-0.7000f,-0.7000f,-0.7000f,
                -0.1000f,-0.1000f,-0.1000f,-0.7000f,-0.7000f,-0.7000f,-0.1000f,-0.1000f,-0.1000f,-0.7000f,-0.7000f,-0.7000f };
            float[] v0 = new float[36] { 
                0.0400f,0.0400f,0.0400f,0.0400f,0.0400f,0.0400f,0.1600f,0.1600f,0.1600f,0.1600f,0.1600f,0.1600f,
                0.0400f,0.0400f,0.0400f,0.0400f,0.0400f,0.0400f,0.1600f,0.1600f,0.1600f,0.1600f,0.1600f,0.1600f,
                0.0400f,0.0400f,0.0400f,0.0400f,0.0400f,0.0400f,0.1600f,0.1600f,0.1600f,0.1600f,0.1600f,0.1600f };
            float[] BNPrice = new float[36] { 
                10.0000f,2.1254f,0.1091f,9.9997f,2.1267f,0.1274f,10.7100f,4.2158f,1.1667f,10.6804f,4.2140f,1.1939f,
                10.1706f,3.4747f,0.7736f,10.1206f,3.4807f,0.8416f,12.1819f,6.4964f,3.0914f,12.1122f,6.4899f,3.1456f,
                10.6478f,4.6473f,1.6832f,10.5637f,4.6636f,1.7874f,13.3142f,8.0083f,4.5454f,13.2172f,7.9998f,4.6201f };

            // Tree settings
            int NT = 100;
            float threshold = 0.01f;

            // Arrays and structures for the American and European prices
            BivariateStruct output = new BivariateStruct();
            float[] EuroPrice = new float[36];
            float[] AmerPrice = new float[36];
            float[] EuroClosedPrice = new float[36];
            float[] CVAmerPrice = new float[36];


            // Create the American and European prices from the tree
            for(int s=0;s<=35;s++)
            {
                output = BuildBivariateTree(S0[s],PutCall,Strike,T[s],rf,NT,kappa,theta,sigma,v0[s],rho[s],threshold);
                EuroPrice[s] = output.EuroPrice;
                AmerPrice[s] = output.AmerPrice;

                // Closed Form European price
                param.v0 = v0[s];
                param.rho = rho[s];
                EuroClosedPrice[s] = HestonPriceGaussLaguerre(param,S0[s],Strike,rf,q,T[s],trap,PutCall,X,W);

                // American price by control variate method
                CVAmerPrice[s] = AmerPrice[s] + (EuroClosedPrice[s] - EuroPrice[s]);
            }

            // Output the results
            Console.WriteLine("---------------------------------------------------------------------------");
            Console.WriteLine("Heston American Option using Believa-Nawalkha tree -- Exhibit 14");
            Console.WriteLine("Using {0:0} time steps",NT);
            Console.WriteLine("---------------------------------------------------------------------------");
            Console.WriteLine(" Spot  rho  sqrt(v0)    T       AmerPutTree  AmerPutCV  TruePrice   EuroPut");
            Console.WriteLine("---------------------------------------------------------------------------");
            for(int s=0;s<=35;s++)
            {
                Console.WriteLine(" {0,3:0} {1,5:F1}  {2,5:F2} {3,10:F4}  {4,10:F5}  {5,10:F5}  {6,10:F5} {7,10:F5}",S0[s],rho[s],Math.Sqrt(v0[s]),T[s],AmerPrice[s],CVAmerPrice[s],BNPrice[s],EuroClosedPrice[s]);
                if ((s+1) % 12 == 0)
                    Console.WriteLine("---------------------------------------------------------------------------");
            }
*/


