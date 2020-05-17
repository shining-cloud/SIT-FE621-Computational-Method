using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Benhamou_Gobet_Miri_Piecewise_DIA_Estimation
{
    class ObjectiveFunction
    {
        // Objective function ===========================================================================
        public double f(double[] param,OFSet ofset)
        {
            BGMPrice BGM = new BGMPrice();
            BisectionAlgo BA = new BisectionAlgo();

            double S = ofset.opset.S;
            double r = ofset.opset.r;
            double q = ofset.opset.q;
            int trap = ofset.opset.trap;
            int ObjFunction = ofset.ObjFunction;
            OPSet opset = ofset.opset;

            double[,] MktIV    = ofset.data.MktIV;
            double[,] MktPrice = ofset.data.MktPrice;
            string[,] PutCall  = ofset.data.PutCall;
            double[] K = ofset.data.K;
            double[] T = ofset.data.T;
            double[] X = ofset.X;
            double[] W = ofset.W;

            int NK = PutCall.GetLength(0);
            int NT = PutCall.GetLength(1);

            // Separate out the parameters from the "param" vector;
            int N = param.Length;
            double kappa = param[0];
            double v0    = param[1];
            double[] THETA = new double[NT];
            double[] SIGMA = new double[NT];
            double[] RHO   = new double[NT];
            for (int k=0; k<=NT-1; k++)
            {
                THETA[k] = param[3*k+2];
                SIGMA[k] = param[3*k+3];
                RHO[k]   = param[3*k+4];
            }

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

            List<double> MatList   = new List<double>();
            List<double> thetaList = new List<double>();
            List<double> sigmaList = new List<double>();
            List<double> rhoList   = new List<double>();

            double[] Mat,theta,sigma,rho;

            if ((kappa<=kappaLB) || (kappa>=kappaUB) || (v0<=v0LB) || (v0>=v0UB) )
                Error = 1e50;
            for (int k=0; k<=NT-1; k++)
                if ((THETA[k]<=thetaLB) || (THETA[k]>=thetaUB) || (SIGMA[k]<=sigmaLB) || (SIGMA[k]>=sigmaUB) || (RHO[k]<=rhoLB) || (RHO[k]>=rhoUB))
                    Error = 1e50;
            {
            for(int t=0;t<NT;t++)
                {
                MatList.Add(T[t]);
                thetaList.Add(THETA[t]);
                sigmaList.Add(SIGMA[t]);
                rhoList.Add(RHO[t]);
                Mat   = MatList.ToArray();
                theta = thetaList.ToArray();
                sigma = sigmaList.ToArray();
                rho   = rhoList.ToArray();
                Mat.Reverse();
                Array.Reverse(theta);
                Array.Reverse(sigma);
                Array.Reverse(rho);

                for(int k=0;k<NK;k++)
                    {
                        ModelPrice[k,t] = BGM.BGMApproxPriceTD(kappa,v0,theta,sigma,rho,opset,K[k],Mat);
                        switch(ObjFunction)
                        {
                            case 1:
                                // MSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2.0);
                                break;
                            case 2:
                                // RMSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2.0) / MktPrice[k,t];
                                break;
                            case 3:
                                // IVMSE Loss Function
                                ModelIV[k,t] = BA.BisecBSIV(opset,K[k],T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                                Error += Math.Pow(ModelIV[k,t] - MktIV[k,t],2);
                                break;
                            case 4:
                                // IVRMSE Christoffersen, Heston, Jacobs proxy
                                double d = (Math.Log(S/K[k]) + (r-q+MktIV[k,t]*MktIV[k,t]/2.0)*T[t])/MktIV[k,t]/Math.Sqrt(T[t]);
                                double NormPDF = Math.Exp(-0.5*d*d)/Math.Sqrt(2.0*pi);
                                Vega = S*NormPDF*Math.Sqrt(T[t]);
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2) / (Vega*Vega);
                                break;
                        }
                    }
                }
            }
            return Error;
        }
    }
}