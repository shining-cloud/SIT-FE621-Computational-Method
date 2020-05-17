using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Chiarella_Ziogas_American_Call
{
    class NewtonMethod
    {
        public double CZNewton(double start,double v0,double v1,double mat,double kappa,double theta,double lambda,double rho,double sigma,
              double K,double r,double q,double[] xs,double[] ws,double[] xt,double[] wt,int Nt,double B0,double B1,int ACNum,double tol,
              double A,double B,double C,double D,string DoubleType)
        {
            // Parameters
            HParam param0;
            param0.kappa = kappa;
            param0.theta = theta;
            param0.sigma = sigma;
            param0.v0 = v0;
            param0.rho = rho;
            param0.lambda = lambda;
            HParam param1;
            param1.kappa = kappa;
            param1.theta = theta;
            param1.sigma = sigma;
            param1.v0 = v1;
            param1.rho = rho;
            param1.lambda = lambda;

            double db = 0.001;
            double diff = 1.1*tol;
            double g0,g,g_,b0,b1,CallA1,CallA0;
            double[] Prices = new double[2];
            double b_new  = 0.0;

            // Set the first value to the starting value
            double b = start;

            // Set the upper time-integration range to the maturity increment
            B = mat;
            CZPrices CZ = new CZPrices();
            while(Math.Abs(diff) > tol)
            {
                if(ACNum == 1)
                {
                    // First set of functions and derivative
                    Prices = CZ.CZAmerCall(Math.Exp(B0 + b*v0),mat,param0,K,r,q,xs,ws,xt,wt,Nt,B0,b,A,B,C,D,DoubleType);
                    CallA1 = Prices[0];
                    g0 = (Math.Log(CallA1 + K) - B0)/v0 - b;

                    Prices = CZ.CZAmerCall(Math.Exp(B0 + (b+db)*v0),mat,param0,K,r,q,xs,ws,xt,wt,Nt,B0,b+db,A,B,C,D,DoubleType);
                    CallA1 = Prices[0];
                    g = (Math.Log(CallA1 + K) - B0)/v0 - (b+db);

                    Prices = CZ.CZAmerCall(Math.Exp(B0 + (b-db)*v0),mat,param0,K,r,q,xs,ws,xt,wt,Nt,B0,b-db,A,B,C,D,DoubleType);
                    CallA1 = Prices[0];
                    g_ = (Math.Log(CallA1 + K) - B0)/v0 - (b-db);
                }
                else
                {
                    // Second set of functions and derivatives
                    Prices = CZ.CZAmerCall(Math.Exp(b + B1*v1),mat,param1,K,r,q,xs,ws,xt,wt,Nt,b,B1,A,B,C,D,DoubleType);
                    CallA0 = Prices[0];
                    g0 = (Math.Log(CallA0 + K) - v1*B1) - b;

                    Prices = CZ.CZAmerCall(Math.Exp(b+db + B1*v1),mat,param1,K,r,q,xs,ws,xt,wt,Nt,b+db,B1,A,B,C,D,DoubleType);
                    CallA0 = Prices[0];
                    g = (Math.Log(CallA0 + K) - v1*B1) - (b+db);

                    Prices = CZ.CZAmerCall(Math.Exp(b-db + B1*v1),mat,param1,K,r,q,xs,ws,xt,wt,Nt,b-db,B1,A,B,C,D,DoubleType);
                    CallA0 = Prices[0];
                    g_ = (Math.Log(CallA0 + K) - v1*B1) - (b-db);
                }
                // The derivative
                double dg = (g - g_)/2.0/db;

                // Newton's method
                b_new = b - g0/dg;
                diff = b_new - b;
                b = b_new;
            }
            return b;
        }
    }
}

