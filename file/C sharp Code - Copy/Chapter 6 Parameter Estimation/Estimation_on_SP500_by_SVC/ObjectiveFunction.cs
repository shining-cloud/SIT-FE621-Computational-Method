using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Estimation_on_SP500_by_SVC
{
    class ObjectiveFunction
    {
        // Objective function ===========================================================================
        public double f(double[] param,OFSet ofset)
        {
            HestonPrice HP = new HestonPrice();
            Bisection BA = new Bisection();

            double S = ofset.opset.S;
            double r = ofset.opset.r;
            double q = ofset.opset.q;
            int trap = ofset.opset.trap;

            double[,] MktIV    = ofset.data.MktIV;
            double[,] MktPrice = ofset.data.MktPrice;
            string[,] PutCall  = ofset.data.PutCall;
            double[] K = ofset.data.K;
            double[] T = ofset.data.T;
            double[] X = ofset.X;
            double[] W = ofset.W;
            string CF = ofset.CF;
            int LossFunction = ofset.LossFunction;

            int NK = PutCall.GetLength(0);
            int NT = PutCall.GetLength(1);
            int NX  = X.Length;

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

            if((param2.kappa<=kappaLB) || (param2.theta<=thetaLB) || (param2.sigma<=sigmaLB) || (param2.v0<=v0LB) || (param2.rho<=rhoLB) || 
               (param2.kappa>=kappaUB) || (param2.theta>=thetaUB) || (param2.sigma>=sigmaUB) || (param2.v0>=v0UB) || (param2.rho>=rhoUB))
                Error = 1e50;
            else
            {
                Complex phi = new Complex(0.0,0.0);
                double phi2 = 0.0;
                Complex i   = new Complex(0.0,1.0);
                Complex[] f2 = new Complex[NX];
                Complex[] f1 = new Complex[NX];
                Complex[] f  = new Complex[NX];
                double[] int1 = new double[NX];
                double[] int2 = new double[NX];
                Complex I1 = new Complex(0.0,0.0);
                Complex I2 = new Complex(0.0,0.0);
                double CallPrice = 0.0;

                for(int t=0;t<NT;t++)
                {
                    for(int j=0;j<NX;j++)
                    {
                        phi = X[j];
                        if(CF == "Heston")
                        {
                            f2[j] = HP.HestonCF(phi,param2,S,r,q,T[t],trap);
                            f1[j] = HP.HestonCF(phi-i,param2,S,r,q,T[t],trap) / (S*Math.Exp((r-q)*T[t]));
                        }
                        else if(CF == "Attari")
                            phi2 = X[j];
                        f[j] = HP.AttariCF(phi2,param2,T[t],S,r,q,trap);
                    }
                    for(int k=0;k<NK;k++)
                    {
                        double L = Math.Log(Math.Exp(-r*T[t])*K[k]/S);
                        for(int j=0;j<NX;j++)
                        {
                            phi = X[j];
                            if(CF == "Heston")
                            {
                                I1 = Complex.Exp(-i*phi*Complex.Log(K[k]))*f1[j] / i / phi;
                                int1[j] = W[j] * I1.Real;
                                I2 = Complex.Exp(-i*phi*Complex.Log(K[k]))*f2[j] / i / phi;
                                int2[j] = W[j] * I2.Real;
                            }
                            else if(CF == "Attari")
                            {
                                phi2 = X[j];
                                double fR = f[j].Real;
                                double fI = f[j].Imaginary;
                                int1[j] = W[j] * ((fR + fI/phi2)*Math.Cos(L*phi2) + (fI - fR/phi2)*Math.Sin(L*phi2)) / (1 + phi2*phi2);
                            }
                        }
                        if(CF == "Heston")
                        {
                            double P1 = 0.5 + 1.0/pi*int1.Sum();
                            double P2 = 0.5 + 1.0/pi*int2.Sum();
                            CallPrice = S*Math.Exp(-q*T[t])*P1 - K[k]*Math.Exp(-r*T[t])*P2;
                        }
                        else if(CF == "Attari")
                            CallPrice = S*Math.Exp(-q*T[t]) - K[k]*Math.Exp(-r*T[t])*(0.5 + 1.0/pi*int1.Sum());

                        if(PutCall[k,t] == "C")
                            ModelPrice[k,t] = CallPrice;
                        else
                            ModelPrice[k,t] = CallPrice - S*Math.Exp(-q*T[t]) + Math.Exp(-r*T[t])*K[k];

                        // Select the objective function
                        switch(LossFunction)
                        {
                            case 1:
                                // MSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2) / Convert.ToDouble(NT*NK); ;
                                break;
                            case 2:
                                // RMSE Loss Function
                                Error += Math.Pow(ModelPrice[k,t] - MktPrice[k,t],2) / MktPrice[k,t] / Convert.ToDouble(NT*NK); ;
                                break;
                            case 3:
                                // IVMSE Loss Function
                                ModelIV[k,t] = BA.BisecBSIV(PutCall[k,t],S,K[k],r,q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                                Error += Math.Pow(ModelIV[k,t] - MktIV[k,t],2) / Convert.ToDouble(NT*NK); ;
                                break;
                            case 4:
                                // IVRMSE Christoffersen, Heston, Jacobs proxy
                                double d = (Math.Log(S/K[k]) + (r-q+MktIV[k,t]*MktIV[k,t]/2.0)*T[t])/MktIV[k,t]/Math.Sqrt(T[t]);
                                double NormPDF = Math.Exp(-0.5*d*d)/Math.Sqrt(2*pi);
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
 
