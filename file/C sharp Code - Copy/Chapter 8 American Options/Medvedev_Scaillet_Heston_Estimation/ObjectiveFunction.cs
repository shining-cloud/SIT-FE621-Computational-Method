using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Medvedev_Scaillet_Heston
{
    class ObjectiveFunction
    {
        // Objective function ===========================================================================
        public double f(double[] param,OFSet ofsettings)
        {
            // Option price settings
            double S = ofsettings.opsettings.S;
            double r = ofsettings.opsettings.r;
            double q = ofsettings.opsettings.q;
            double T = ofsettings.opsettings.T;
            int trap = ofsettings.opsettings.trap;

            // Market data
            double[] MktIV = ofsettings.data.MktIV;
            int NK = MktIV.Length;

            // MS settings
            double dt   = ofsettings.mssettings.dt;
            double a    = ofsettings.mssettings.a;
            double b    = ofsettings.mssettings.b;
            double tol  = ofsettings.mssettings.tol;
            int MaxIter = ofsettings.mssettings.MaxIter;

            // Optimizatin settings
            double[] lb = ofsettings.lb;
            double[] ub = ofsettings.ub;
            double[] X  = ofsettings.X;
            double[] W  = ofsettings.W;

            HParam param2 = new HParam();
            param2.kappa = param[0];
            param2.theta = param[1];
            param2.sigma = param[2];
            param2.v0    = param[3];
            param2.rho   = param[4];
            param2.lambda = 0.0;

            // Settings for the Bisection algorithm
            a = 0.01;
            b = 3.0;
            double B = 3.0;

            // Initialize the model price and model implied vol vectors, and the objective function value
            double[] ModelPrice = new double[NK];
            double[] ModelIV    = new double[NK];
            double[] Error      = new double[NK];
            double SumError = 0.0;

            // Parameter bounds
            double kappaLB = lb[0]; double kappaUB = ub[0];
            double thetaLB = lb[1]; double thetaUB = ub[1];
            double sigmaLB = lb[2]; double sigmaUB = ub[2];
            double v0LB    = lb[3]; double v0UB    = ub[3];
            double rhoLB   = lb[4]; double rhoUB   = ub[4];

            // Classes
            MSPrices MS = new MSPrices();
            BisectionImpliedVol BA = new BisectionImpliedVol();

            // Penalty for inadmissible parameter values
            if((param2.kappa<=kappaLB) || (param2.theta<=thetaLB) || (param2.sigma<=sigmaLB) || (param2.v0<=v0LB) || (param2.rho<=rhoLB) || 
               (param2.kappa>=kappaUB) || (param2.theta>=thetaUB) || (param2.sigma>=sigmaUB) || (param2.v0>=v0UB) || (param2.rho>=rhoUB))
                SumError = 1.0e50;

            // Penalty for inadmissible implied vol
            else
            {
                Console.WriteLine("----------------------------------------");
                Console.WriteLine("ModelPrice ImpliedVol   MktVol    Error");
                Console.WriteLine("----------------------------------------");
                for(int k=0;k<=NK-1;k++)
                {
                    double Strike = ofsettings.data.K[k];
                    ofsettings.opsettings.K = Strike;
                    double[] output = new double[6];
                    output = MS.MSPriceHeston(param2,ofsettings.opsettings,ofsettings.mssettings);
                    ModelPrice[k] = Math.Max(0.01,output[2]);
                    ModelIV[k] = BA.BisectionMSIV(S,Strike,r,q,T,a,b,ModelPrice[k],tol,MaxIter,B,dt);
                    if(ModelIV[k] == -1.0)
                        Error[k] = 1.0e50;
                    else
                        Error[k] = Math.Pow(ModelIV[k] - MktIV[k],2.0);
                    Console.WriteLine("{0,7:F4} {1,10:F4} {2,10:F4} {3,10:F6}",ModelPrice[k],ModelIV[k],MktIV[k],Error[k]);
                    SumError += Error[k];
                }
            }
            return SumError;
        }
    }
}
