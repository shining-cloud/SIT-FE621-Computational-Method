using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Heston
{
    class MSPrices
    {
        public double MSPriceBS(double S,double K,double T,double sigma,double r,double q,int MaxIter,double Tol,double b,double dt)
        {
            MSsetBS mssettings = new MSsetBS();
            mssettings.theta = Math.Log(K/S)/sigma/Math.Sqrt(T);
            mssettings.K = K;
            mssettings.sigma = sigma;
            mssettings.r = r;
            mssettings.q = q;
            mssettings.T = T;
            Tol = 1e-5;
            MaxIter = 50;
            double a = Math.Max(1.25,0.95*mssettings.theta);
            b = 3.0;

            // Golden section search method to find y
            GoldenSearch GS = new GoldenSearch();
            double y = GS.GoldenSearchBS(a,b,Tol,MaxIter,mssettings);

            // Bisection algorithm to find y
            //double y = BisectionBS(mssettings,a,b,Tol,MaxIter,dt);


            if(y<mssettings.theta)
                y = mssettings.theta;
            MSExpansionBS MS = new MSExpansionBS();
            return MS.MSPutBS(y,mssettings);
        }
        public double[] MSPriceHeston(HParam param,OpSet opset,MSSet msset)
        {
            // param  = Heston parameters
            // opset  = option settings
            // msset  = Medvedev-Scaillet expansion settings

            // Extract quantities from the option settings
            double K = opset.K;
            double S = opset.S;
            double v0 = param.v0;
            double T = opset.T;
            double r = opset.r;
            double q = opset.q;

            // Extract quantities from the MS structure
            int N = msset.N;
            double a = msset.a;
            double b = msset.b;
            double A = msset.A;
            double B = msset.B;
            double dt  = msset.dt;
            double tol = msset.tol;
            int MaxIter  = msset.MaxIter;
            int NumTerms = msset.NumTerms;
            double yinf  = msset.yinf;
            int method   = msset.method;

            // Closed-form European put
            NewtonCotes NC = new NewtonCotes();
            double EuroPutClosed = NC.HestonPriceNewtonCotes(param,opset,method,A,B,N);
            
            // Moneyness
            double theta = Math.Log(K/S)/Math.Sqrt(v0)/Math.Sqrt(T);

            // Find the barrier level
            // Bisection method to find y
            // a = Math.Max(1.25,0.95*theta);
            // b = 3.0;
            // double y = BisectionMS(a,b,theta,K,param,r,q,T,NumTerms,tol,MaxIter,dt);

            // Golden section search method to find y
            GoldenSearch GS = new GoldenSearch();
            a = 1.5;
            b = 2.5;
            double y = GS.GoldenSearchMS(a,b,tol,MaxIter,K,param,theta,r,q,T,NumTerms);
            if(y<theta)
                y = theta;

            // Find the early exercise premium
            MSExpansionHeston MS = new MSExpansionHeston();
            double AmerPutMS = MS.MSPutHeston(y,theta,K,param,r,q,T,NumTerms);
            double EuroPutMS = MS.MSPutHeston(yinf,theta,K,param,r,q,T,NumTerms);
            double EEP = AmerPutMS - EuroPutMS;
            
            // Control variate American put
            double AmerPut = EuroPutClosed + EEP;
            
            // Output the results
            double[] output = new double[6];
            output[0] = EuroPutClosed;
            output[1] = AmerPutMS;
            output[2] = AmerPut;      // Medvedev and Scaillet recommend this one.
            output[3] = EEP;
            output[4] = theta;
            output[5] = y;
            return output;
        }
    }
}

