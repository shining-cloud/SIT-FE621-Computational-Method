using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Atiya_Wall_MLE
{
    class Likelihood
    {
        public double f(double[] param,OFSet ofset)
        {
            // Name the Heston parameters
            double kappa = param[0];
            double theta = param[1];
            double sigma = param[2];
            double v0    = param[3];
            double rho   = param[4];
             
            // Likelihood function settings
            double[] x = ofset.x;
            double r = ofset.r;
            double q = ofset.q;
            double dt = ofset.dt;
            int Lmethod = ofset.method;

            // Atiya and Wall parameterization
            double alpha = kappa*theta;
            double beta  = kappa;

            // Number of log-stock prices
            int T = x.Length;

            // Drift term
            double mu = r - q;

            // Equation (17)
            double betap = 1.0 - beta*dt;

            // Equation (18) - denominator of d(t)
            double D = 2.0*Math.PI*sigma*Math.Sqrt(1.0-rho*rho)*dt;

            // Equation (14)
            double a = (betap*betap + rho*sigma*betap*dt + sigma*sigma*dt*dt/4.0) / (2.0*sigma*sigma*(1.0-rho*rho)*dt);

            // Variance and likelihood at time t = 0
            double[] v = new Double[T];
            double[] L = new Double[T];
            v[0] = v0;
            if(Lmethod==1)
                L[0] = Math.Exp(-v[0]);    // Construct the Likelihood
            else if(Lmethod==2)
                L[0] = -v[0];         // Construct the log-likelihood

            // Construction the likelihood for time t = 1 through t = T
            double dx,B,C,bt,x1,x2,E;
            for(int t=0;t<=T-2;t++)
            {
                // Stock price increment
                dx  = x[t+1] - x[t];
                // Equations (31) and (32)
                B = -alpha*dt - rho*sigma*(dx-mu*dt);
                C = alpha*alpha*dt*dt + 2.0*rho*sigma*alpha*dt*(dx-mu*dt) + sigma*sigma*Math.Pow(dx-mu*dt,2.0) - 2.0*v[t]*v[t]*a*sigma*sigma*(1.0-rho*rho)*dt;
                // Equation (30) to update the variance
                if(B*B - C > 0.0)
                    v[t+1] = Math.Sqrt(B*B - C) - B;
                else
                {
                    // If v[t+1] is negative, use the approximation Equation (33)
                    bt = (Math.Pow(v[t]-alpha*dt,2.0) - 2.0*rho*sigma*(v[t]-alpha*dt)*(dx-mu*dt) + sigma*sigma*Math.Pow(dx-mu*dt,2.0))  / (2.0*sigma*sigma*(1.0-rho*rho)*dt);
                    if(bt/a > 0.0)
                        v[t+1] = Math.Sqrt(bt/a);
                    else
                        // If v[t+1] is still negative, take the previous value
                        v[t+1] = v[t];
                }
                // Equation (15) and (16)
                bt = (Math.Pow(v[t+1]-alpha*dt,2.0) - 2.0*rho*sigma*(v[t+1]-alpha*dt)*(dx-mu*dt) + sigma*sigma*Math.Pow(dx-mu*dt,2.0))  / (2.0*sigma*sigma*(1.0-rho*rho)*dt);
                x1 = ((2.0*betap+rho*sigma*dt)*(v[t+1]-alpha*dt) - (2.0*rho*sigma*betap+sigma*sigma*dt)*(dx-mu*dt))   / (2.0*sigma*sigma*(1.0-rho*rho)*dt);
                x2 = -2.0*Math.Sqrt(a*bt);
                // Compbined exponent for Equation (34)
                E = Math.Exp(x1 + x2) / D;
                if(Lmethod==1)
                    // Equation (34) for the likelihood L[t+1]
                    L[t+1] = Math.Pow(a*bt,-0.25) * E * L[t];
                else if(Lmethod==2)
                    // Alternatively, use the log-likelihood, log of Equation (34)
                    L[t+1] = -0.25*Math.Log(a*bt) + x1 + x2 - Math.Log(D) + L[t];
            }
            // Negative likelihood is the last term.
            // Since we maximize the likelihood, we minimize the negative likelihood.
            return -L[T-1];
        }
    }
}


