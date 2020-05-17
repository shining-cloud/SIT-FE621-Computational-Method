using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Quadratic_Exponential
{
    class QESimulation
    {
        // Price by simulation
        public double QEPrice(HParam param,OpSet settings,double gamma1,double gamma2,int NT,int NS,int MC,double phic,string PutCall)
        {
            RandomNumbers RN = new RandomNumbers();
            double[] STe = QESim(param,settings,gamma1,gamma2,NT,NS,MC,phic);
            double[] Price = new double[NS];
            for(int s=0;s<=NS-1;s++)
            {
                if(settings.PutCall == "C")
                    Price[s] = Math.Max(STe[s] - settings.K,0.0);
                else if(settings.PutCall == "P")
                    Price[s] = Math.Max(settings.K - STe[s],0.0);
            }
            return Math.Exp(-settings.r*settings.T) * RN.Mean(Price);
        }

        // Simulation of stock price paths and variance paths using Euler or Milstein schemes
        public double[] QESim(HParam param,OpSet settings,double gamma1,double gamma2,int NT,int NS,int MC,double phic)
        {
            // Heston parameters
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            // option settings
            double r = settings.r;
            double q = settings.q;

            // Time increment
            double dt = settings.T/Convert.ToDouble(NT);

            // Initialize the variance and stock processes
            double[,] V = new double[NT,NS];
            double[,] S = new double[NT,NS];

            // Flags for negative variances
            int F = 0;

            // Starting values for the variance and stock processes
            for(int s=0;s<=NS-1;s++)
            {
                S[0,s] = settings.S;       // Spot price 
                V[0,s] = v0;       // Heston v0 initial variance 
            }

            // Parameters for the QE scheme
            double E = Math.Exp(-kappa*dt);
            double K1 = (kappa*rho/sigma - 0.5)*gamma1*dt - rho/sigma;
            double K2 = (kappa*rho/sigma - 0.5)*gamma2*dt + rho/sigma;
            double K3 = gamma1*dt*(1.0 - rho*rho);
            double K4 = gamma2*dt*(1.0 - rho*rho);
            double A = K2 + K4/2.0;
            double m,s2,phi,Uv,b,a,K0,p,beta,M;
            double phiinv = 0.0;

            // Generate the stock and volatility paths
            RandomNumbers RN = new RandomNumbers();
            InverseNorm IN = new InverseNorm();
            double Zv,Zs;
            for(int s=0;s<=NS-1;s++)
            {
                for(int t=1;t<=NT-1;t++)
                {
                    // Generate two dependent N(0,1) variables with correlation rho
                    Zv = RN.RandomNorm();
                    Zs = rho*Zv + Math.Sqrt(1.0-rho*rho)*RN.RandomNorm();

                    // QE Agorithm
                    m = theta + (V[t-1,s] - theta)*E;
                    s2 = V[t-1,s]*sigma*sigma*E/kappa*(1-E) + theta*sigma*sigma/2.0/kappa*(1.0-E)*(1.0-E);
                    phi = s2/m/m;
                    Uv = RN.RandomNum(0.0,1.0);
                    if(phi <= phic)
                    {
                        b = Math.Sqrt(2.0/phi - 1.0 + Math.Sqrt(2.0/phi*(2.0/phi-1.0)));
                        a = m/(1.0+b*b);
                        Zv = IN.normICDF(Uv);
                        V[t,s] = a*(b + Zv)*(b + Zv);
                        K0 = -kappa*rho*theta*dt/sigma;

                        // Martingale correction: Define new K0 if possible
                        if((MC==1) & A<(1.0/(2.0*a)))
                        {
                            M = Math.Exp(A*b*b*a/(1.0-2.0*A*a)) / Math.Sqrt(1.0-2.0*A*a);
                            K0 = -Math.Log(M) - (K1 + 0.5*K3)*V[t-1,s];
                        }

                        // Stock price
                        S[t,s] = S[t-1,s]*Math.Exp((r-q)*dt + K0 + K1*V[t-1,s] + K2*V[t,s] + Math.Sqrt(K3*V[t-1,s] + K4*V[t,s])*Zs);
                    }
                    else
                    {
                        p = (phi-1.0)/(phi+1.0);
                        beta = (1.0-p)/m;
                        if((0.0<=Uv) & (Uv<=p))
                            phiinv = 0.0;
                        else if((p<Uv) & (Uv<=1.0))
                            phiinv = 1.0/beta*Math.Log((1.0-p)/(1.0-Uv));
                        V[t,s] = phiinv;
                        K0 = -kappa*rho*theta*dt/sigma;

                        // Martingale correction: Define new K0 if possible
                        if(MC==1 & A<beta)
                        {
                            M = p + beta*(1.0-p)/(beta-A);
                            K0 = -Math.Log(M) - (K1 + 0.5*K3)*V[t-1,s];
                        }

                        // Stock price
                        S[t,s] = S[t-1,s]*Math.Exp((r-q)*dt + K0 + K1*V[t-1,s] + K2*V[t,s] + Math.Sqrt(K3*V[t-1,s] + K4*V[t,s])*Zs);
                    }
                }
            }
            // Return the vector or terminal stock prices
            double[] output = new double[NS];
            for(int s=0;s<=NS-1;s++)
                output[s] = S[NT-1,s];

            return output;
        }
    }
}

