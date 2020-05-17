using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Mikhailov_and_Nogel_Estimation_DJIA
{
    class ObjectiveFunction
    {
        // Objective function ===========================================================================
        public double f(double[] param,double tau,double[] tau0,double[,] param0,MktData data,double S,double[] K,double r,double q,int trap,
                        double[] X,double[] W,int LossFunction,double[] lb,double[] ub)
        {
            HestonPriceTD HPTD = new HestonPriceTD();
            BisectionAlgo BA = new BisectionAlgo();
            
            double[] MktIV    = data.MktIV;
            double[] MktPrice = data.MktPrice;
            string[] PutCall  = data.PutCall;

            int NK = PutCall.Length;
            int NT = 1;

            HParam param2 = new HParam();
            param2.kappa = param[0];
            param2.theta = param[1];
            param2.sigma = param[2];
            param2.v0    = param[3];
            param2.rho   = param[4];

            // Settings for the Bisection algorithm
            double a = 0.001;
            double b = 3.0;
            double Tol = 1e5;
            int MaxIter = 10000;

            // Initialize the model price and model implied vol vectors, and the objective function value
            double[] ModelPrice = new double[NK];
            double[] ModelIV    = new double[NK];
            double Vega = 0.0;
            double Error = 0.0;
            double pi = Math.PI;

            double kappaLB = lb[0]; double kappaUB = ub[0];
            double thetaLB = lb[1]; double thetaUB = ub[1];
            double sigmaLB = lb[2]; double sigmaUB = ub[2];
            double v0LB    = lb[3]; double v0UB    = ub[3];
            double rhoLB   = lb[4]; double rhoUB   = ub[4];

            if((param2.kappa<=kappaLB) || (param2.theta<=thetaLB) || (param2.sigma<=sigmaLB) || (param2.v0<=v0LB) || (param2.rho<=rhoLB) || 
               (param2.kappa>=kappaUB) || (param2.theta>=thetaUB) || (param2.sigma>=sigmaUB) || (param2.v0>=v0UB) || (param2.rho>=rhoUB))
                Error = 1e50;
            else
            {
                for(int k=0;k<=NK-1;k++)
                {
                    ModelPrice[k] = HPTD.MNPriceGaussLaguerre(param2,param0,tau,tau0,S,K[k],r,q,PutCall[k],trap,X,W);
                    switch(LossFunction)
                    {
                        case 1:
                            // MSE Loss Function
                            Error += Math.Pow(ModelPrice[k] - MktPrice[k],2);
                            break;
                        case 2:
                            // RMSE Loss Function
                            Error += Math.Pow(ModelPrice[k] - MktPrice[k],2) / MktPrice[k];
                            break;
                        case 3:
                            // IVMSE Loss Function
                            ModelIV[k] = BA.BisecBSIV(PutCall[k],S,K[k],r,q,tau,a,b,ModelPrice[k],Tol,MaxIter);
                            Error += Math.Pow(ModelIV[k] - MktIV[k],2);
                            break;
                        case 4:
                            // IVRMSE Christoffersen, Heston, Jacobs proxy
                            double d = (Math.Log(S/K[k]) + (r-q+MktIV[k]*MktIV[k]/2.0)*tau)/MktIV[k]/Math.Sqrt(tau);
                            double NormPDF = Math.Exp(-0.5*d*d)/Math.Sqrt(2.0*pi);
                            Vega = S*NormPDF*Math.Sqrt(tau);
                            Error += Math.Pow(ModelPrice[k] - MktPrice[k],2) / (Vega*Vega);
                            break;
                    }
                }
            }
            return Error;  //       / Convert.ToDouble(NT*NK);
        }
    }
}
