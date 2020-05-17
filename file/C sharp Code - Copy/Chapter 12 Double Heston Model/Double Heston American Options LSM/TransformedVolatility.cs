using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Double_Heston_American_Options_LSM
{
    class TVSimulation
    {
        public DHSim DHTransVolSim(string scheme,DHParam param,double S0,double Strike,double Mat,double r,double q,int T,int N,string PutCall)
        {
            // Heston parameters
            double kappa1 = param.kappa1;
            double theta1 = param.theta1;
            double sigma1 = param.sigma1;
            double v01    = param.v01;
            double rho1   = param.rho1;
            double kappa2 = param.kappa2;
            double theta2 = param.theta2;
            double sigma2 = param.sigma2;
            double v02    = param.v02;
            double rho2   = param.rho2;

            // Time increment
            double dt = Mat/T;

            // Required quantities
            double K01 = -rho1*kappa1*theta1*dt/sigma1;
            double K11 = dt/2.0*(kappa1*rho1/sigma1 - 0.5) - rho1/sigma1;
            double K21 = dt/2.0*(kappa1*rho1/sigma1 - 0.5) + rho1/sigma1;
            double K31 = dt/2.0*(1.0-rho1*rho1);

            double K02 = -rho2*kappa2*theta2*dt/sigma2;
            double K12 = dt/2.0*(kappa2*rho2/sigma2 - 0.5) - rho2/sigma2;
            double K22 = dt/2.0*(kappa2*rho2/sigma2 - 0.5) + rho2/sigma2;
            double K32 = dt/2.0*(1.0-rho2*rho2);

            // Initialize the volatility and stock processes
            double[,] w1 = new double[T,N];
            double[,] w2 = new double[T,N];
            double[,] S  = new double[T,N];

            // Starting values for the variance and stock processes
            for(int k=0;k<=N-1;k++)
            {
                S[0,k]  = S0;                   // Spot price 
                w1[0,k] = Math.Sqrt(v01);       // Heston initial volatility
                w2[0,k] = Math.Sqrt(v02);
            }
            // Generate the stock and volatility paths
            RandomNumber RN = new RandomNumber();
            double Zv1,m11,m12,beta1,thetav1;
            double Zv2,m21,m22,beta2,thetav2;
            double B1,B2,logS;
            for(int k=0;k<=N-1;k++)
                for(int t=1;t<=T-1;t++)
                {
                    Zv1 = RN.RandomNorm();
                    Zv2 = RN.RandomNorm();
                    if(scheme == "ZhuEuler")
                    {
                        // Transformed variance scheme
                        w1[t,k] = w1[t-1,k] + 0.5*kappa1*((theta1-sigma1*sigma1/4.0/kappa1)/w1[t-1,k] - w1[t-1,k])*dt + 0.5*sigma1*Math.Sqrt(dt)*Zv1;
                        w2[t,k] = w2[t-1,k] + 0.5*kappa2*((theta2-sigma2*sigma2/4.0/kappa2)/w2[t-1,k] - w2[t-1,k])*dt + 0.5*sigma2*Math.Sqrt(dt)*Zv2;
                    }
                    else if(scheme == "ZhuTV")
                    {
                        // Zhu (2010) process for the transformed volatility
                        m11 = theta1 + (w1[t-1,k]*w1[t-1,k] - theta1)*Math.Exp(-kappa1*dt);
                        m21 = theta2 + (w2[t-1,k]*w2[t-1,k] - theta2)*Math.Exp(-kappa2*dt);
                        m12 = sigma1*sigma1/4.0/kappa1*(1.0 - Math.Exp(-kappa1*dt));
                        m22 = sigma2*sigma2/4.0/kappa2*(1.0 - Math.Exp(-kappa2*dt));
                        beta1 = Math.Sqrt(Math.Max(0,m11-m12));
                        beta2 = Math.Sqrt(Math.Max(0,m21-m22));
                        thetav1 = (beta1 - w1[t-1,k]*Math.Exp(-kappa1*dt/2.0))/(1.0 - Math.Exp(-kappa1*dt/2.0));
                        thetav2 = (beta2 - w2[t-1,k]*Math.Exp(-kappa2*dt/2.0))/(1.0 - Math.Exp(-kappa2*dt/2.0));
                        w1[t,k] = w1[t-1,k] + 0.5*kappa1*(thetav1 - w1[t-1,k])*dt + 0.5*sigma1*Math.Sqrt(dt)*Zv1;
                        w2[t,k] = w2[t-1,k] + 0.5*kappa2*(thetav2 - w2[t-1,k])*dt + 0.5*sigma2*Math.Sqrt(dt)*Zv2;
                    }
                    // Predictor-Corrector for the stock price
                    B1 = RN.RandomNorm();
                    B2 = RN.RandomNorm();
                    logS = Math.Log(Math.Exp(-r*Convert.ToDouble(t)*dt)*S[t-1,k])
			             + K01 + K11*w1[t-1,k]*w1[t-1,k] + K21*w1[t,k]*w1[t,k] + Math.Sqrt(K31*(w1[t,k]*w1[t,k] + w1[t-1,k]*w1[t-1,k]))*B1
			             + K02 + K12*w2[t-1,k]*w2[t-1,k] + K22*w2[t,k]*w2[t,k] + Math.Sqrt(K32*(w2[t,k]*w2[t,k] + w2[t-1,k]*w2[t-1,k]))*B2;
                    S[t,k] = Math.Exp(logS)*Math.Exp(r*Convert.ToDouble(t+1)*dt);
                }
            // Terminal stock prices
            double[] ST = new double[N];
            for(int k=0;k<=N-1;k++)
                ST[k] = S[T-1,k];

            // Payoff vectors
            double[] Payoff = new double[N];
            for(int k=0;k<=N-1;k++)
            {
                if(PutCall == "C")
                    Payoff[k] = Math.Max(ST[k] - Strike,0.0);
                else if(PutCall == "P")
                    Payoff[k] = Math.Max(Strike - ST[k],0.0);
            }
            // Simulated price
            Regression RE = new Regression();
            double SimPrice = Math.Exp(-r*Mat)*RE.VMean(Payoff);

            // Output the results
            DHSim output = new DHSim();
            output.S = S;
            output.EuroPrice = SimPrice;

            return output;
        }
    }
}

