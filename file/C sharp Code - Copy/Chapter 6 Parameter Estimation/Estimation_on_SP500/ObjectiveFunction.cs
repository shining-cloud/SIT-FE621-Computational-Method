using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Estimation_on_SP500
{
    class ObjectiveFunction
    {
        // Objective function ===========================================================================
        public double f(double[] param,OFSet ofset)
        {
            Bisection BA = new Bisection();
            HestonPrice HP = new HestonPrice();
            double S = ofset.opsettings.S;
            double r = ofset.opsettings.r;
            double q = ofset.opsettings.q;
            int trap = ofset.opsettings.trap;

            double[,] MktIV    = ofset.data.MktIV;
            double[,] MktPrice = ofset.data.MktPrice;
            string[,] PutCall  = ofset.data.PutCall;
            double[] K = ofset.data.K;
            double[] T = ofset.data.T;

            int NK = PutCall.GetLength(0);
            int NT = PutCall.GetLength(1);

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
            double[,] ModelPrice = new double[NK,NT];
            double[,] ModelIV    = new double[NK,NT];
            double Vega = 0.0;
            double Error = 0.0;
            double pi = Math.PI;

            double[] lb = ofset.lb;
            double[] ub = ofset.ub;

            double kappaLB = lb[0]; double kappaUB = ub[0];
            double thetaLB = lb[1]; double thetaUB = ub[1];
            double sigmaLB = lb[2]; double sigmaUB = ub[2];
            double v0LB    = lb[3]; double v0UB    = ub[3];
            double rhoLB   = lb[4]; double rhoUB   = ub[4];

            int LossFunction = ofset.LossFunction;
            double[] X = ofset.X;
            double[] W = ofset.W;

            // Penalty for inadmissible parameter values
            if((param2.kappa<=kappaLB) || (param2.theta<=thetaLB) || (param2.sigma<=sigmaLB) || (param2.v0<=v0LB) || (param2.rho<=rhoLB) || 
               (param2.kappa>=kappaUB) || (param2.theta>=thetaUB) || (param2.sigma>=sigmaUB) || (param2.v0>=v0UB) || (param2.rho>=rhoUB))
                Error = 1e50;
            else
            {
                for(int k=0;k<NK;k++)
                {
                    for(int t=0;t<NT;t++)
                    {
                        ModelPrice[k,t] = HP.HestonPriceGaussLaguerre(param2,S,K[k],r,q,T[t],trap,PutCall[k,t],X,W);
                        switch(LossFunction)
                        {
                            case 1:
                                // MSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2) / Convert.ToDouble(NT*NK);
                                break;
                            case 2:
                                // RMSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2) / MktPrice[k,t] / Convert.ToDouble(NT*NK);
                                break;
                            case 3:
                                // IVMSE Loss Function
                                ModelIV[k,t] = BA.BisecBSIV(PutCall[k,t],S,K[k],r,q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                                Error += Math.Pow(ModelIV[k,t] - MktIV[k,t],2) / Convert.ToDouble(NT*NK); ;
                                break;
                            case 4:
                                // IVRMSE Christoffersen, Heston, Jacobs proxy
                                double d = (Math.Log(S/K[k]) + (r-q+MktIV[k,t]*MktIV[k,t]/2.0)*T[t])/MktIV[k,t]/Math.Sqrt(T[t]);
                                double NormPDF = Math.Exp(-0.5*d*d)/Math.Sqrt(2.0*pi);
                                Vega = S*NormPDF*Math.Sqrt(T[t]);
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2) / Vega / Vega / Convert.ToDouble(NT*NK);
                                break;
                        }
                    }
                }
            }
            return Error;
        }
    }
}