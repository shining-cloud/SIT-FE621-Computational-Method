using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Elices_TimeDependent_DIA_Estimation
{
    class ObjectiveFunction
    {
        // Objective function for Elices prices ===========================================================================
        public double f(double[] param,OFSetE ofsettingsE,OFSetH ofsettingsH,string choice)
        {
            ElicesAlgo EA = new ElicesAlgo();
            HestonPrice HP = new HestonPrice();

            double Error = 0.0;
            if(choice=="Elices")
            {
                double S = ofsettingsE.opsettings.S;
                double r = ofsettingsE.opsettings.r;
                double q = ofsettingsE.opsettings.q;
                int trap = ofsettingsE.opsettings.trap;

                double[] MktIV    = ofsettingsE.data.MktIV;
                double[] MktPrice = ofsettingsE.data.MktPrice;
                string PutCall    = ofsettingsE.data.PutCall;
                double[] K = ofsettingsE.data.K;
                int NK = K.Length;
                int t = ofsettingsE.t;

                double[][] paramfixed = ofsettingsE.paramfixed;
                double[] T = ofsettingsE.T;
                double v0  = ofsettingsE.v0;
                int NT = T.Length;

                double[] lb = ofsettingsE.lb;
                double[] ub = ofsettingsE.ub;

                int LossFunction = ofsettingsE.LossFunction;
                double[] X = ofsettingsE.X;
                double[] W = ofsettingsE.W;

                // Initialize the model price and model implied vol vectors, and the objective function value
                double[] ModelPrice = new double[NK];
                double[] ModelIV    = new double[NK];
                double Vega = 0.0;
                double pi = Math.PI;

                double kappaLB = lb[0]; double kappaUB = ub[0];
                double thetaLB = lb[1]; double thetaUB = ub[1];
                double sigmaLB = lb[2]; double sigmaUB = ub[2];
                double rhoLB   = lb[3]; double rhoUB   = ub[3];

                // Penalty for inadmissible parameter values
                if((param[0]<=kappaLB) || (param[1]<=thetaLB) || (param[2]<=sigmaLB) || (param[3]<=rhoLB) || 
                   (param[0]>=kappaUB) || (param[1]>=thetaUB) || (param[2]>=sigmaUB) || (param[3]>=rhoUB))
                    Error = 1.0e50;
                else
                {
                    for(int k=0;k<NK-1;k++)
                    {
                        if(t==0)
                            ModelPrice[k] = EA.ElicesPrice(PutCall,S,K[k],T,r,q,param,v0,trap,X,W);
                        else
                            ModelPrice[k] = EA.ElicesPrice(PutCall,S,K[k],T,r,q,param,v0,trap,X,W,paramfixed);
                        switch(LossFunction)
                        {
                            case 1:
                                // MSE Loss Function
                                Error += Math.Pow(ModelPrice[k] - MktPrice[k],2.0);
                                break;
                            case 2:
                                // RMSE Loss Function
                                Error += Math.Pow(ModelPrice[k] - MktPrice[k],2.0) / MktPrice[k];
                                break;
                            case 3:
                                // IVRMSE Christoffersen, Heston, Jacobs proxy
                                double mat = T[NT-1];
                                double d = (Math.Log(S/K[k]) + (r+MktIV[k]*MktIV[k]/2.0)*mat)/MktIV[k]/Math.Sqrt(mat);
                                double NormPDF = Math.Exp(-0.5*d*d)/Math.Sqrt(2.0*pi);
                                Vega = S*NormPDF*Math.Sqrt(mat);
                                Error += Math.Pow(ModelPrice[k] - MktPrice[k],2.0) / (Vega*Vega);
                                break;
                        }
                    }
                    Error /= Convert.ToDouble(NK);
                }
            }
            else if(choice == "Heston")
            // Objective function for Heston prices ===========================================================================
            {
                double[,] MktPrice = ofsettingsH.MktPrice;
                double[,] MktIV    = ofsettingsH.MktIV;
                string PutCall    = ofsettingsH.PutCall;
                double[] K = ofsettingsH.K;
                int NK = K.Length;

                double[] T = ofsettingsH.T;
                int NT = T.Length;

                double[] lb = ofsettingsH.lb;
                double[] ub = ofsettingsH.ub;

                int LossFunction = ofsettingsH.LossFunction;
                double[] X = ofsettingsH.X;
                double[] W = ofsettingsH.W;

                // Initialize the model price and model implied vol vectors, and the objective function value
                double[,] ModelPrice = new double[NK,NT];
                double[,] ModelIV    = new double[NK,NT];
                double Vega = 0.0;
                double pi = Math.PI;

                double kappaLB = lb[0]; double kappaUB = ub[0];
                double thetaLB = lb[1]; double thetaUB = ub[1];
                double sigmaLB = lb[2]; double sigmaUB = ub[2];
                double rhoLB   = lb[3]; double rhoUB   = ub[3];
                double v0LB    = lb[4]; double v0UB    = ub[4];

                HParam param2 = new HParam();
                param2.kappa = param[0];
                param2.theta = param[1];
                param2.sigma = param[2];
                param2.rho   = param[3];
                param2.v0    = param[4];

                // Penalty for inadmissible parameter values
                if((param[0]<=kappaLB) || (param[1]<=thetaLB) || (param[2]<=sigmaLB) || (param[3]<=rhoLB) || (param[4]<=v0LB) ||
                   (param[0]>=kappaUB) || (param[1]>=thetaUB) || (param[2]>=sigmaUB) || (param[3]>=rhoUB) || (param[4]>=v0UB))
                    Error = 1.0e50;
                else
                {
                    for(int t=0;t<=NT-1;t++)
                    {
                        for(int k=0;k<=NK-1;k++)
                        {
                            ModelPrice[k,t] = HP.HestonPriceGaussLaguerre(param2,ofsettingsH.opsettings,K[k],T[t],PutCall,X,W);
                            switch(LossFunction)
                            {
                                case 1:
                                    // MSE Loss Function
                                    Error += Math.Pow(ModelPrice[k,t] - MktPrice[t,k],2.0);
                                    break;
                                case 2:
                                    // RMSE Loss Function
                                    Error += Math.Pow(ModelPrice[k,t] - MktPrice[t,k],2.0) / MktPrice[k,t];
                                    break;
                                case 3:
                                    // IVRMSE Christoffersen, Heston, Jacobs proxy
                                    double S = ofsettingsH.opsettings.S;
                                    double r = ofsettingsH.opsettings.r;
                                    double q = ofsettingsH.opsettings.q;
                                    double mat = T[NT-1];
                                    double d = (Math.Log(S/K[k]) + (r+MktIV[t,k]*MktIV[t,k]/2.0)*mat)/MktIV[k,t]/Math.Sqrt(mat);
                                    double NormPDF = Math.Exp(-0.5*d*d)/Math.Sqrt(2.0*pi);
                                    Vega = S*NormPDF*Math.Sqrt(mat);
                                    Error += Math.Pow(ModelPrice[k,t] - MktPrice[t,k],2.0) / (Vega*Vega);
                                    break;
                            }
                        }
                    }
                    Error /= Convert.ToDouble(NK*NT);
                }
            }
            return Error;
        }
    }
}
