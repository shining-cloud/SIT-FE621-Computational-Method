using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Double_Heston_Simulation
{
    class EulerAlfonsiSimulation
    {
        public double DHEulerAlfonsiSim(string scheme,DHParam param,double S0,double Strike,double Mat,double r,double q,int T,int N,string PutCall)
        {
            RandomNumber RN = new RandomNumber();

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

            HParam param1;
            param1.kappa = kappa1;
            param1.theta = theta1;
            param1.sigma = sigma1;
            param1.v0    = v01;
            param1.rho   = rho1;

            HParam param2;
            param2.kappa = kappa2;
            param2.theta = theta2;
            param2.sigma = sigma2;
            param2.v0    = v02;
            param2.rho   = rho2;

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
            double G1,G2,B1,B2,logS;

            // Generate the stock and volatility paths
            for(int k=0;k<=N-1;k++)
                for(int t=1;t<=T-1;t++)
                {
                    if(scheme == "Euler")
                    {
                        // Generate two in dependent N(0,1) variables
                        G1 = RN.RandomNorm();
                        G2 = RN.RandomNorm();
                        // Euler discretization with full truncation for the variances
                        V1[t,k] = V1[t-1,k] + kappa1*(theta1-V1[t-1,k])*dt + sigma1*Math.Sqrt(V1[t-1,k]*dt)*G1;
                        V2[t,k] = V2[t-1,k] + kappa2*(theta2-V2[t-1,k])*dt + sigma2*Math.Sqrt(V2[t-1,k]*dt)*G2;
                        V1[t,k] = Math.Max(0,V1[t,k]);
                        V2[t,k] = Math.Max(0,V2[t,k]);
                    }
                    else if(scheme  == "Alfonsi")
                    {
                        // Alfonsi discretization
                        V1[t,k] = AlfonsiV(param1,V1[t-1,k],dt);
                        V2[t,k] = AlfonsiV(param2,V2[t-1,k],dt);
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
            double SimPrice = Math.Exp(-r*Mat)*VMean(Payoff);
            return SimPrice;
        }

        // Alfonsi function
        public double AlfonsiV(HParam param,double vt,double dt)
        {
            RandomNumber RN = new RandomNumber();
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;

            double phi = (1.0-Math.Exp(-kappa*dt/2.0))/kappa;
            double S = (sigma*sigma/4.0 - theta*kappa);
            double E = Math.Exp(kappa*dt/2.0);

            double K2 = 0.0;
            double U = 0.0;
            double Y = 0.0;
            double newV = 0.0;
            if(sigma*sigma > 4.0*kappa*theta)
                K2 = E*(S*phi + Math.Pow((Math.Sqrt(E*S*phi) + sigma/2.0*Math.Sqrt(3.0*dt)),2));
            else
                K2 = 0;
            if(vt >= K2)
            {
                U = RN.RandomNum(0.0,1.0);
                if(U <= 1.0/6.0)
                    Y = Math.Sqrt(3.0);
                else if(U <= 1.0/3.0)
                    Y = -Math.Sqrt(3.0);
                else
                    Y = 0.0;
                phi = (1.0-Math.Exp(-kappa*dt/2.0))/kappa;
                S = (theta*kappa - sigma*sigma/4.0);
                E = Math.Exp(-kappa*dt/2.0);
                newV = E*Math.Pow((Math.Sqrt(S*phi + E*vt) + sigma/2.0*Math.Sqrt(dt)*Y),2) + S*phi;
            }
            else
            {
                double[] u = CIRmoments(param,vt,dt);
                double u1 = u[0];
                double u2 = u[1];
                double Pi = 0.5 - 0.5*Math.Sqrt(1 - u1*u1/u2);
                U = RN.RandomNum(0.0,1.0);
                if(U <= Pi)
                    newV = u1/2.0/Pi;
                else if(U > Pi)
                    newV = u1/2.0/(1.0-Pi);
            }
            return newV;
        }
        // CIR moments
        public double[] CIRmoments(HParam param,double Vs,double dt)
        {
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;

            // E[vt | vs];
            double e = theta + (Vs - theta)*Math.Exp(-kappa*dt);

            // Var[vt | vs]
            double v = Vs*sigma*sigma*Math.Exp(-kappa*dt)/kappa*(1.0-Math.Exp(-kappa*dt))
              + theta*sigma*sigma/2.0/kappa*Math.Pow((1.0-Math.Exp(-kappa*dt)),2);

            // E[vt^2 | vs]
            double e2 = v + e*e;

            double[] output = new double[2];
            output[0] = e;
            output[1] = e2;
            return output;
        }

        // Mean of a vector
        public double VMean(double[] X)
        {
            double XBar = 0.0;
            int N = X.Length;
            for(int i=0;i<=N-1;i++)
                XBar += X[i];
            return XBar / Convert.ToDouble(N);
        }
    }
}



