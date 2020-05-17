using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Alfonsi
{
    class AlfonsiSim
    {
        public double AlfonsiPrice(HParam param,double S0,double Strike,double Mat,double r,double q,int T,int N,string PutCall)
        {
            RandomNumber RN = new RandomNumber();

            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;

            // Time increment
            double dt = Mat/T;

            // Required quantities
            double K0 = -rho*kappa*theta*dt/sigma;
            double K1 = dt/2.0*(kappa*rho/sigma - 0.5) - rho/sigma;
            double K2 = dt/2.0*(kappa*rho/sigma - 0.5) + rho/sigma;
            double K3 = dt/2.0*(1.0-rho*rho);

            // Initialize the variance and stock processes
            double[,] V = new double[T,N];
            double[,] S = new double[T,N];

            // Starting values for the variance and stock processes
            for(int k=0;k<=N-1;k++)
            {
                S[0,k] = S0;       // Spot price 
                V[0,k] = v0;       // Heston v0 initial variance 
            }
            double B,logS;

            // Generate the stock and volatility paths
            for(int k=0;k<=N-1;k++)
                for(int t=1;t<=T-1;t++)
                {
                    // Alfonsi discretization
                    V[t,k] = AlfonsiV(param,V[t-1,k],dt);
                    // Predictor-Corrector for the stock price
                    B = RN.RandomNorm();
                    logS = Math.Log(Math.Exp(-r*Convert.ToDouble(t)*dt)*S[t-1,k])
			             + K0 + K1*V[t-1,k] + K2*V[t,k] + Math.Sqrt(K3*(V[t,k]+V[t-1,k]))*B;
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

            // K2 parameter
            if(sigma*sigma > 4.0*kappa*theta)
                K2 = E*(S*phi + Math.Pow((Math.Sqrt(E*S*phi) + sigma/2.0*Math.Sqrt(3.0*dt)),2.0));
            else
                K2 = 0.0;

            // Update the variance
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
                newV = E*Math.Pow((Math.Sqrt(S*phi + E*vt) + sigma/2.0*Math.Sqrt(dt)*Y),2.0) + S*phi;
            }
            else
            {
                double[] u = CIRmoments(param,vt,dt);
                double u1 = u[0];
                double u2 = u[1];
                double Pi = 0.5 - 0.5*Math.Sqrt(1.0 - u1*u1/u2);
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
              + theta*sigma*sigma/2.0/kappa*Math.Pow((1.0-Math.Exp(-kappa*dt)),2.0);

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

