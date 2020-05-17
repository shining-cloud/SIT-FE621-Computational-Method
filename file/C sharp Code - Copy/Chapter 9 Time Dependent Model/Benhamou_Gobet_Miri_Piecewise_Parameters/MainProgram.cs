using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Benhamou_Gobet_Miri_Piecewise_Parameters
{
    class BGM_Piecewise
    {
        static void Main(string[] args)
        {
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

            // Reproduces ATM prices in Table 8 of Benhamou, Gobet, and Miri
            //"Time Dependent Heston Model"
            // SIAM Journal on Financial Mathematics, Vol. 1, (2010)

            // Exact prices from the BGM (2010) paper
            double[] TruePW = new double[8] { 3.93,5.53,7.85,11.23,13.92,18.37,22.15,27.17 };

            // Spot price, risk free rate, dividend yield
            OpSet settings = new OpSet();
            settings.S  = 100.0;
            settings.r = 0.0;
            settings.q  = 0.0;
            settings.PutCall = "P";
            settings.trap = 1;       // Little trap formulation
            double K = 100.0;

            // Heston piecewise constant parameters.
            double kappa = 3.0;
            double v0 = 0.04;
            double[] T = new double[40];
            double[] THETA = new double[40];
            double[] SIGMA = new double[40];
            double[] RHO   = new double[40];

            // Construct the piecewise constant parameters
            for(int t=0;t<=39;t++)
            {
                T[t]     = (Convert.ToDouble(t)+1)/4.0;
                THETA[t] =  0.04 + Convert.ToDouble(t)*0.05/100.0;
                SIGMA[t] =  0.30 + Convert.ToDouble(t)*0.50/100.0;
                RHO[t]   = -0.20 + Convert.ToDouble(t)*0.35/100.0;
            }

            // Compute maturity intervals for the Mikhailov-Nogel model
            double[] tau = new double[40];
            tau[0] = T[0];
            for(int t=1;t<=39;t++)
                tau[t] = T[t] - T[t-1];

            // Obtain the Approximate Piecewise and Closed Form Piecewise put prices using the constructed parameters
            List<double> MatList    = new List<double>();
            List<List<double>> MNparam0 = new List<List<double>>();
            List<double> tau0       = new List<double>();
            List<double> thetaList  = new List<double>();
            List<double> sigmaList  = new List<double>();
            List<double> rhoList    = new List<double>();

            HParam param = new HParam();
            BGMPrice BGM = new BGMPrice();
            HestonPrice HP = new HestonPrice();
            HestonPriceTD HPTD = new HestonPriceTD();

            double[] ApproxPW = new double[40];
            double[] ClosedPW = new double[40];
            double[] Mat,theta,sigma,rho;
            for (int t=0; t<=39; t++)
            {
                // Stack the maturities and parameters, with the oldest ones at the top, and the newest ones the bottom
                MatList.Add(T[t]);
                Mat = MatList.ToArray();
                thetaList.Add(THETA[t]);
                sigmaList.Add(SIGMA[t]);
                rhoList.Add(RHO[t]);
                theta = thetaList.ToArray();
                sigma = sigmaList.ToArray();
                rho   = rhoList.ToArray();
                // Calculate the approximate formula with the piecewise constant parameters
                ApproxPW[t] = BGM.BGMApproxPriceTD(kappa,v0,theta,sigma,rho,settings,K,Mat);
                 if(t==0)
                {
                    // First iteration for the closed-form price with the PW constant parameters
                    param.kappa = kappa;
                    param.theta = theta[0];
                    param.sigma = sigma[0];
                    param.v0    = v0;
                    param.rho   = rho[0];
                    ClosedPW[t] = HP.HestonPriceGaussLaguerre(param,settings,K,T[0],X,W);
                }
                else if(t==1)
                {
                    // Second iteration for the closed form price with the PW constant parameters
                    double[] OldMat = new double[1] { 0.25 };
                    double[,] param0 = new double[1,5];
                    param0[0,0] = kappa;
                    param0[0,1] = theta[0]; param0[0,2] = sigma[0]; param0[0,3] = v0; param0[0,4] = rho[0];
                    ClosedPW[t] = HPTD.HestonPriceGaussLaguerreTD(param,param0,tau[t],OldMat,settings,K,X,W);
                }
                else if(t > 1)
                {
                    // Remaing iteration for the closed form price with PW constant parameters
                    // Current parameter values
                    param.theta = theta[t];
                    param.sigma = sigma[t];
                    param.rho   = rho[t];
                    // Stack the past maturities, oldest on the bottm, newest on top
                    double[] OldMat = new double[t];
                    for(int k=0;k<=t-1;k++)
                        OldMat[k] = tau[t-k-1];
                    // Fill in the past parameter values, with the oldest on the bottom and newest on top
                    double[,] param0 = new double[t,5];
                    for(int k=0;k<=t-1;k++)
                    {
                        param0[k,0] = kappa;
                        param0[k,1] = theta[t-k-1];
                        param0[k,2] = sigma[t-k-1];
                        param0[k,3] = v0;
                        param0[k,4] = rho[t-k-1];
                    }
                    ClosedPW[t] = HPTD.HestonPriceGaussLaguerreTD(param,param0,tau[t],OldMat,settings,K,X,W);
                }
            }

            // Select only the puts at the maturities
            double[] ClosedPW2 = new double[8];
            double[] ApproxPW2 = new double[8];
            double[] Mats = new double[8] {0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0};

            for (int t=0; t<=39; t++)
                for(int k=0;k<=7;k++)
                    if(T[t] == Mats[k])
                    {
                        ClosedPW2[k] = ClosedPW[t];
                        ApproxPW2[k] = ApproxPW[t];
                    }

            // The averaged parameter values
            // Note: v0(4) has been changed
            double[] v0new    = new double[] { 0.04,  0.0397,  0.0328,  0.0464,  0.05624, 0.2858,  0.8492,  0.1454};
            double[] thetanew = new double[] { 0.04,  0.0404,  0.0438,  0.0402,  0.0404,  0.0268,  0.0059,  0.0457};
            double[] sigmanew = new double[] { 0.30,  0.3012,  0.3089,  0.3112,  0.3210,  0.3363,  0.3541,  0.3998};
            double[] rhonew   = new double[] {-0.20, -0.1993, -0.1972, -0.1895, -0.1820, -0.1652, -0.1480, -0.1232};

           // The closed form prices using the averaged parameter volumes
            HParam newparam = new HParam();
            double[] ClosedAvg = new double[8];
            for(int t=0;t<=7;t++)
            {
                newparam.kappa = kappa;
                newparam.theta = thetanew[t];
                newparam.sigma = sigmanew[t];
                newparam.v0    = v0new[t];
                newparam.rho   = rhonew[t];
                ClosedAvg[t] = HP.HestonPriceGaussLaguerre(newparam,settings,K,Mats[t],X,W);
            }

            // Output the results at the maturities
            Console.WriteLine("----------------------------------------------------------------------------");
            Console.WriteLine("          Original Heston | Mihailov-Nogel |  BGM Matlab    | BGM True Price");
            Console.WriteLine("Maturity    Closed Form   |   Piecewise    |   Piecewise    |   Piecewise   ");
            Console.WriteLine("----------------------------------------------------------------------------");
            for(int t=0;t<=7;t++)
                Console.WriteLine("{0,5:0.00} {1,14:F2} {2,16:F2} {3,16:F2} {4,16:F2}",Mats[t],ClosedAvg[t],ClosedPW2[t],ApproxPW2[t],TruePW[t]);
            Console.WriteLine("----------------------------------------------------------------------------");
        }
    }
}

