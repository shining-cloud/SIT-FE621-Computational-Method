using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Double_Heston_Estimation
{
    class ObjectiveFunction
    {
        // Objective function ===========================================================================
        public double DHObjFun(double[] param,OpSet settings,MktData data,double[] X,double[] W,int objectivefun,double[] lb,double[] ub)
        {
            HestonPriceDH HPDH = new HestonPriceDH();
            BisectionAlgo BA = new BisectionAlgo();

            double S = settings.S;
            double r = settings.r;
            double q = settings.q;
            int trap = settings.trap;

            double[,] MktIV    = data.MktIV;
            double[,] MktPrice = data.MktPrice;
            string[,] PutCall  = data.PutCall;
            double[] K = data.K;
            double[] T = data.T;

            int NK = PutCall.GetLength(0);
            int NT = PutCall.GetLength(1);

            DHParam param2 = new DHParam();
            param2.kappa1  = param[0];
            param2.theta1  = param[1];
            param2.sigma1  = param[2];
            param2.v01     = param[3];
            param2.rho1    = param[4];
            param2.kappa2  = param[5];
            param2.theta2  = param[6];
            param2.sigma2  = param[7];
            param2.v02     = param[8];
            param2.rho2    = param[9];
            
            // Settings for the Bisection algorithm
            double a = 0.001;
            double b = 3.0;
            double Tol = 1e5;
            int MaxIter = 10000;

            // Initialize the model price and model implied vol vectors, and the objective function value
            double[,] ModelPrice = new double[NK,NT];
            double[,] ModelIV    = new double[NK,NT];
            double BSVega = 0.0;
            double Error = 0.0;
            double pi = Math.PI;

            double kappaLB = lb[0]; double kappaUB = ub[0];
            double thetaLB = lb[1]; double thetaUB = ub[1];
            double sigmaLB = lb[2]; double sigmaUB = ub[2];
            double    v0LB = lb[3]; double    v0UB = ub[3];
            double   rhoLB = lb[4]; double   rhoUB = ub[4]; 

            if((param2.kappa1<=kappaLB) || (param2.theta1<=thetaLB) || (param2.sigma1<=sigmaLB) || (param2.v01<=v0LB) || (param2.rho1<=rhoLB) || 
               (param2.kappa1>=kappaUB) || (param2.theta1>=thetaUB) || (param2.sigma1>=sigmaUB) || (param2.v01>=v0UB) || (param2.rho1>=rhoUB) ||
               (param2.kappa2<=kappaLB) || (param2.theta2<=thetaLB) || (param2.sigma2<=sigmaLB) || (param2.v02<=v0LB) || (param2.rho2<=rhoLB) || 
               (param2.kappa2>=kappaUB) || (param2.theta2>=thetaUB) || (param2.sigma2>=sigmaUB) || (param2.v02>=v0UB) || (param2.rho2>=rhoUB))
                Error = 1e50;
            else
            {
                for(int k=0;k<NK;k++)
                {
                    for(int t=0;t<NT;t++)
                    {
                        settings.K = K[k];
                        settings.T = T[t];
                        ModelPrice[k,t] = HPDH.DoubleHestonPriceGaussLaguerre(param2,settings,X,W);
                        switch(objectivefun)
                        {
                            case 1:
                                // MSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2.0);
                                break;
                            case 2:
                                // RMSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2.0) / MktPrice[k,t] ;
                                break;
                            case 3:
                                // IVMSE Loss Function
                                ModelIV[k,t] = BA.BisecBSIV(PutCall[k,t],S,K[k],r,q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                                Error += Math.Pow(ModelIV[k,t] - MktIV[k,t],2.0) ;
                                break;
                            case 4:
                                // IVRMSE Christoffersen, Heston, Jacobs proxy
                                double d = (Math.Log(S/K[k]) + (r-q+MktIV[k,t]*MktIV[k,t]/2.0)*T[t])/MktIV[k,t]/Math.Sqrt(T[t]);
                                double NormPDF = Math.Exp(-0.5*d*d)/Math.Sqrt(2.0*pi);
                                BSVega = S*NormPDF*Math.Sqrt(T[t]);
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2.0) / (BSVega*BSVega);
                                break;
                        }
                    }
                }
            }
            return Error;
        }
    }
}
