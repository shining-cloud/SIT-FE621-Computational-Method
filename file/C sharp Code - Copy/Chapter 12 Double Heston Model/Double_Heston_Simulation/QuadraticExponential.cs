using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Double_Heston_Simulation
{
    class QESimulation
    {
        public double DHQuadExpSim(DHParam param,double S0,double Strike,double Mat,double r,double q,int T,int N,string PutCall)
        {
            RandomNumber RN = new RandomNumber();
            NormInverse NI = new NormInverse();

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

            // Initialize the variance and stock processes
            double[,] V1 = new double[T,N];
            double[,] V2 = new double[T,N];
            double[,] S  = new double[T,N];

            // Starting values for the variance and stock processes
            for(int k=0;k<=N-1;k++)
            {
                S[0,k]  = S0;       // Spot price 
                V1[0,k] = v01;       // Heston v0 initial variance 
                V2[0,k] = v02;       // Heston v0 initial variance 
            }

            // Generate the stock and volatility paths
            double m1,s1,phi1,p1,U1,b1,a1,Zv1;
            double m2,s2,phi2,p2,U2,b2,a2,Zv2;
            double beta1,beta2,B1,B2,logS;

            for(int k=0;k<=N-1;k++)
                for(int t=1;t<=T-1;t++)
                {
                    m1 = theta1 + (V1[t-1,k] - theta1)*Math.Exp(-kappa1*dt);
                    s1 = V1[t-1,k]*sigma1*sigma1*Math.Exp(-kappa1*dt)/kappa1*(1.0-Math.Exp(-kappa1*dt))
                       + theta1*sigma1*sigma1/2.0/kappa1*Math.Pow(1.0-Math.Exp(-kappa1*dt),2.0);
                    phi1 = s1/(m1*m1);
                    p1 = (phi1-1.0)/(phi1+1.0);
                    U1 = RN.RandomNum(0.0,1.0);
                    if(phi1 < 0.5)
                    {
                        b1 = Math.Sqrt(2.0/phi1 - 1.0 + Math.Sqrt(2.0/phi1*(2.0/phi1-1.0)));
                        a1 = m1/(1.0 + b1*b1);
                        Zv1 = NI.normICDF(U1);
                        V1[t,k] = a1*(b1+Zv1)*(b1+Zv1);
                    }
                    else if(phi1 >= 0.5)
                    {
                        if(U1 <= p1)
                            V1[t,k] = 0.0;
                        else if(U1 > p1)
                        {
                            beta1 = (1.0-p1)/m1;
                            V1[t,k] = Math.Log((1.0-p1)/(1-U1))/beta1;
                        }
                    }
                    m2 = theta2 + (V2[t-1,k] - theta2)*Math.Exp(-kappa2*dt);
                    s2 = V2[t-1,k]*sigma2*sigma2*Math.Exp(-kappa2*dt)/kappa2*(1.0-Math.Exp(-kappa2*dt))
                       + theta2*sigma2*sigma2/2.0/kappa2*Math.Pow(1.0-Math.Exp(-kappa2*dt),2.0);
                    phi2 = s2/(m2*m2);
                    p2 = (phi2-1.0)/(phi2+1.0);
                    U2 = RN.RandomNum(0.0,1.0);
                    if(phi2 < 0.5)
                    {
                        b2 = Math.Sqrt(2.0/phi2 - 1.0 + Math.Sqrt(2.0/phi2*(2.0/phi2-1.0)));
                        a2 = m2/(1.0 + b2*b2);
                        Zv2 = NI.normICDF(U2);
                        V2[t,k] = a2*(b2+Zv2)*(b2+Zv2);
                    }
                    else if(phi2 >= 0.5)
                    {
                        if(U2 <= p2)
                            V2[t,k] = 0.0;
                        else if(U2 > p2)
                        {
                            beta2 = (1.0-p2)/m2;
                            V2[t,k] = Math.Log((1.0-p2)/(1.0-U2))/beta2;
                        }
                    }
                    // Predictor-Corrector for the stock price
                    B1 = RN.RandomNorm();
                    B2 = RN.RandomNorm();
                    logS = Math.Log(Math.Exp(-r*Convert.ToDouble(t)*dt)*S[t-1,k])
			             + K01 + K11*V1[t-1,k] + K21*V1[t,k] + Math.Sqrt(K31*(V1[t,k]+V1[t-1,k]))*B1
			             + K02 + K12*V2[t-1,k] + K22*V2[t,k] + Math.Sqrt(K32*(V2[t,k]+V2[t-1,k]))*B2;
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
            EulerAlfonsiSimulation EA = new EulerAlfonsiSimulation();
            double SimPrice = Math.Exp(-r*Mat)*EA.VMean(Payoff);
            return SimPrice;
        }
    }
}



