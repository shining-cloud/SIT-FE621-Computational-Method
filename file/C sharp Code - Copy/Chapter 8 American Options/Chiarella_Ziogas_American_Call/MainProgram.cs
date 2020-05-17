using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Chiarella_Ziogas_American_Call
{
    class CZAmericanCall
    {
        static void Main(string[] args)
        {
            // 32-point Gauss-Laguerre Abscissas and weights
            double[] xs = new double[32];
            double[] ws = new double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    xs[k] = double.Parse(bits[0]);
                    ws[k] = double.Parse(bits[1]);
                }
            // 32-point Gauss-Legendre Abscissas and weights
            double[] xt = new double[32];
            double[] wt = new double[32];
            using(TextReader reader = File.OpenText("../../GaussLegendre32.txt"))
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    xt[k] = double.Parse(bits[0]);
                    wt[k] = double.Parse(bits[1]);
                }

            // Option settings
            double S0 = 100.0;
            double K  = 100.0;
            double tau = 0.25;
            double r = 0.01;
            double q = 0.12;
            int trap = 1;

            // Option settings into the structure
            OPSet opset = new OPSet();
            opset.S = S0;
            opset.K = K;
            opset.r = r;
            opset.q = q;

            // Heston parameters
            HParam param = new HParam();
            param.kappa  = 4.0;
            param.theta  = 0.09;
            param.sigma  = 0.1;
            param.v0     = 0.04;
            param.rho    = 0.0;
            param.lambda = 0.0;

            // Time integration limits
            double a = 1e-10;
            double b = tau;
            
            // Price integration limits
            double c = 1e-10;
            double d = 150;

            // Number of points for double trapezoidal rule
            int Nt = 50;

            // Choose double integration "GLe" or "Trapz"
            string DoubleType = "GLe";

            // Closed form European price
            HestonPrices HP = new HestonPrices();
            double EuroClosed = HP.HestonPriceGaussLaguerre(param,opset,"C",tau,trap,xs,ws);

            // Starting values for v0[n], v1[n], b0[n], and b1[n]
            findBfunctions FB = new findBfunctions();
            double Evt = FB.EV(param.v0,param.theta,param.kappa,tau,tau);
            double V00 = Evt + param.sigma/Math.Abs(param.kappa)*Math.Sqrt(param.kappa*param.theta/2.0);
            double V10 = Evt - param.sigma/Math.Abs(param.kappa)*Math.Sqrt(param.kappa*param.theta/2.0);
            double b00 = Math.Log(K);
            double b10 = 0.0;

            // Tolerance and maximum number of steps
            double tol0 = 0.005;
            double tol1 = 0.005;
            int Ntau = 25;
            double Ntol = 1e-8;

            // Vectors B0[n] and B1[n] from the Chiarella algorithm
            BVec bvector = FB.findB(tau,param,K,r,q,V00,V10,b00,b10,xs,ws,xt,wt,Nt,Ntau,tol0,tol1,Ntol,a,b,c,d,DoubleType);
            double[] B0 = bvector.B0;
            double[] B1 = bvector.B1;
            int Nb = B0.Length;

            // Last values of the vectors go into the American Call price
            double b0 = B0[Nb-1];
            double b1 = B1[Nb-1];

            // Compute the call prices with b0 and b1 returned from the optimization
            CZPrices CZ = new CZPrices();
            double[] CallPrices = new double[2];
            CallPrices = CZ.CZAmerCall(S0,tau,param,K,r,q,xs,ws,xt,wt,Nt,b0,b1,a,b,c,d,DoubleType);
            double AmerCZ = CallPrices[0];
            double EuroCZ = CallPrices[1];

            // Prices from the Explicit Method, wUsing a Non-Uniform Grid
            // Minimum and maximum values for the Stock Price, Volatility, and Maturity
            double Smin = 0.0; double Smax = 2.0*K;
            double Vmin = 0.0; double Vmax = 0.5;
            double Tmin = 0.0; double Tmax = tau;

            // Number of grid points for the stock, volatility, and maturity
            int nS = 89;        // Stock price
            int nV = 49;        // Volatility
            int nT = 5000;      // Maturity

            // The maturity time increment and grid
            double dt = (Tmax-Tmin)/Convert.ToDouble(nT);
            double[] T = new double[nT+1];
            for(int i=0;i<=nT;i++)
                T[i] = Convert.ToDouble(i)*dt;

            // The stock price grid
            Interpolation IP = new Interpolation();
            double cc = K/5.0;
            double dz = 1.0/nS*(IP.aSinh((Smax-K)/cc) - IP.aSinh(-K/cc));
            double[] z = new double[nS+1];
            double[] S = new double[nS+1];
            for(int i=0;i<=nS;i++)
            {
                z[i] = IP.aSinh(-K/cc) + Convert.ToDouble(i)*dz;
                S[i] = K + cc*Math.Sinh(z[i]);
            }
            S[0] = 0;

            // The volatility grid
            double dd = Vmax/500.0;
            double dn = IP.aSinh(Vmax/dd)/nV;
            double[] n = new double[nV+1];
            double[] V = new double[nV+1];
            for(int j=0;j<=nV;j++)
            {
                n[j] = Convert.ToDouble(j)*dn;
                V[j] = dd*Math.Sinh(n[j]);
            }
            // Solve the PDE using the Explicit Method
            //double[,] AmerU = HestonExplicitPDENonUniformGrid(param,K,r,q,S,V,T,"C","A");
            //double[,] EuroU = HestonExplicitPDENonUniformGrid(param,K,r,q,S,V,T,"C","E");

            // Obtain the American and European price by 2-D interpolation
            //double AmerPDE = interp2(V,S,AmerU,param.v0,S0);
            //double EuroPDE = interp2(V,S,EuroU,param.v0,S0);

            // Obtain the control variate price from the PDE
            double AmerPDE = 3.730064;
            double EuroPDE = 3.508536;
            double AmerCV = EuroClosed + (AmerPDE - EuroPDE);

            // LSM prices
            //AmerLSM = 3.729049;
            //EuroLSM = 3.501443;
            //AmerCV  = 3.733422;

            // Output the results
            Console.WriteLine("------------------------------------------------------------");
            Console.WriteLine("      b0        b1");
            Console.WriteLine(" {0,10:F5} {1,10:F5}",b0,b1);
            Console.WriteLine("-----------------------------");
            Console.WriteLine("                European   American   AmericanCV");
            Console.WriteLine("-------------------------------------------------");
            Console.WriteLine("Closed        {0,10:F6}                    ",EuroClosed);
            Console.WriteLine("Explicit      {0,10:F6} {1,10:F6} {2,10:F6}",EuroPDE,AmerPDE,AmerCV);
            Console.WriteLine("Chiarella     {0,10:F6} {1,10:F6}          ",EuroCZ, AmerCZ);
            Console.WriteLine("-------------------------------------------------");
        }
    }
}
