using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_Heston
{
    class MSPriceHeston
    {
        // Medvedev-Scaillet price using expansion
        public double[] MSPrice(HParam param,OpSet opset,MSSet msset)
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
            BisectionAlgo BA = new BisectionAlgo();
            MSPut MS = new MSPut();
            double EuroPutClosed = NC.HestonPriceNewtonCotes(param,opset,method,A,B,N);

            // Moneyness
            double theta = Math.Log(K/S)/Math.Sqrt(v0)/Math.Sqrt(T);

            // Find the barrier level
            double y = BA.Bisection(a,b,theta,K,param,r,q,T,NumTerms,tol,MaxIter,dt);
            if(y<theta)
                y = theta;

            // Find the early exercise premium
            double AmerPutMS = MS.MSPutHeston(y,theta,K,param,r,q,T,NumTerms);
            double EuroPutMS = MS.MSPutHeston(yinf,theta,K,param,r,q,T,NumTerms);
            double EEP = AmerPutMS - EuroPutMS;

            // Control variate American put
            double AmerPut = EuroPutClosed + EEP;

            // Output the results
            double[] output = new double[6];
            output[0] = EuroPutClosed;
            output[1] = AmerPutMS;
            output[2] = AmerPut;
            output[3] = EEP;
            output[4] = theta;
            output[5] = y;
            return output;
        }
    }
}

