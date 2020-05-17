using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Elices_TimeDependent_DIA_Estimation
{
    class ElicesAlgo
    {
        // Time dependent characteristic function with fixed parameters
        public Complex ElicesCF(Complex phi1,double[] param,double v0,double[] tau,
                                double S,double r,double q,int trap,double[][] paramfixed)
        {
            Complex i = new Complex(0.0,1.0);  // Imaginary unit
            double x0 = Math.Log(S);

            // phi2(N) coefficient vector
            int N = tau.Length;
            Complex[] phi2 = new Complex[N];
            Complex phi20 = 0.0;
            phi2[N-1] = 0.0;
            if(N>=2)
                phi2[N-2] = -i*C(phi1,phi2[N-1],param,tau[N-1],S,trap);
            if(N>=3)
                for(int t=N-3;t>=0;t--)
                {
                    double[] paramfixed2 = new double[4];
                    paramfixed2[0] = paramfixed[t+1][0];
                    paramfixed2[1] = paramfixed[t+1][1];
                    paramfixed2[2] = paramfixed[t+1][2];
                    paramfixed2[3] = paramfixed[t+1][3];
                    phi2[t] = -i*C(phi1,phi2[t+1],paramfixed2,tau[t+1],S,trap);
                }
            if(N>=2)
            {
                double[] paramfixed2 = new double[4];
                paramfixed2[0] = paramfixed[0][0];
                paramfixed2[1] = paramfixed[0][1];
                paramfixed2[2] = paramfixed[0][2];
                paramfixed2[3] = paramfixed[0][3];

                phi20 = -i*C(phi1,phi2[0],paramfixed2,tau[0],S,trap);
            }
            else
                phi20 = -i*C(phi1,phi2[0],param,tau[0],S,trap);

            // A coefficients.  
            Complex zero = new Complex(0.0,0.0);
            Complex[] Ah = new Complex[N];
            Ah[N-1] = A(phi1,zero,param,tau[N-1],r,q,trap);  // Current params
            if(N>=2)
                for(int t=N-2;t>=0;t--)
                {
                    double[] paramfixed2 = new double[4];
                    paramfixed2[0] = paramfixed[t][0];
                    paramfixed2[1] = paramfixed[t][1];
                    paramfixed2[2] = paramfixed[t][2];
                    paramfixed2[3] = paramfixed[t][3];
                    Ah[t] = A(phi1,phi2[t],paramfixed2,tau[t],r,q,trap);
                }
            // C coefficient
            Complex Ch = new Complex();
            if(N>=2)
            {
                double[] paramfixed2 = new double[4];
                paramfixed2[0] = paramfixed[0][0];
                paramfixed2[1] = paramfixed[0][1];
                paramfixed2[2] = paramfixed[0][2];
                paramfixed2[3] = paramfixed[0][3];
                Ch = C(phi1,phi20,paramfixed2,tau[0],S,trap);
            }
            else
                Ch = C(phi1,phi20,param,tau[0],S,trap);

            Complex AhSum = 0.0;
            for(int j=0;j<=Ah.Length-1;j++)
                AhSum += Ah[j];
            
            return Complex.Exp(AhSum + i*phi1*x0 + Ch*v0);
        }
        // Time dependent characteristic function with NO fixed parameters
        public Complex ElicesCF(Complex phi1,double[] param,double v0,double[] tau,
                                double S,double r,double q,int trap)
        {
            Complex i = new Complex(0.0,1.0);  // Imaginary unit
            double x0 = Math.Log(S);

            // phi2(N) coefficient vector
            int N = tau.Length;
            double Mat = tau[0];
            Complex phi2  = new Complex(0.0,0.0);
            Complex phi20 = -i*C(phi1,phi2,param,Mat,S,trap);

            // A coefficients.  
            Complex zero = new Complex(0.0,0.0);
            Complex Ah = A(phi1,zero,param,Mat,r,q,trap);  // Current params

            // C coefficient
            Complex Ch = C(phi1,phi20,param,Mat,S,trap);
            return Complex.Exp(Ah + i*phi1*x0 + Ch*v0);
        }
        // Time dependent "A" function
        public Complex A(Complex phi1,Complex phi2,double[] param,double tau,double r,double q,int trap)
        {
            Complex i = new Complex(0.0,1.0);                   // Imaginary unit
            double kappa = param[0];
            double theta = param[1];
            double sigma = param[2];
            double rho   = param[3];
            double a = kappa * theta;
            double lambda = 0.0;

            Complex u = -0.5;
            Complex b = kappa + lambda;
            Complex d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi1 - b,2.0) - sigma*sigma*(2.0*u*i*phi1 - phi1*phi1));
            Complex num = (b - rho*sigma*i*phi1 + d - sigma*sigma*i*phi2);
            Complex den = (b - rho*sigma*i*phi1 - d - sigma*sigma*i*phi2);
            if(den==0)
                den = 0.01;
            Complex g = num/den;
            Complex c = 1.0/g;
            Complex G=0.0,AA=0.0;
            if(trap==1)
            {
                // Little Trap formulation in Kahl (2008)
                G = (c*Complex.Exp(-d*tau)-1.0)/(c-1.0);
                AA = (r-q)*i*phi1*tau + a/sigma/sigma*((b - rho*sigma*i*phi1 - d)*tau - 2.0*Complex.Log(G));
            }
            else if(trap==0)
            {
                // Original Heston formulation.
                G = (1.0-g*Complex.Exp(d*tau))/(1.0-g);
                AA = (r-q)*i*phi1*tau + a/sigma/sigma*((b - rho*sigma*i*phi1 + d)*tau - 2.0*Complex.Log(G));
            }
            return AA;
        }
        // Time dependent C function
        public Complex C(Complex phi1,Complex phi2,double[] param,double tau,double S,int trap)
        {
            Complex i = new Complex(0.0,1.0);                   // Imaginary unit
            double kappa = param[0];
            double theta = param[1];
            double sigma = param[2];
            double rho   = param[3];

            double a = kappa * theta;
            double lambda = 0.0;

            Complex u = -0.5;
            Complex b = kappa + lambda;
            Complex d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi1 - b,2.0) - sigma*sigma*(2.0*u*i*phi1 - phi1*phi1));
            Complex g = (b - rho*sigma*i*phi1 + d - sigma*sigma*i*phi2) / (b - rho*sigma*i*phi1 - d - sigma*sigma*i*phi2);
            Complex c = 1.0/g;
            Complex G=0.0,CC=0.0;
            if(trap==1)
            {
                // Little Trap formulation in Kahl (2008)
                G = (c*Complex.Exp(-d*tau)-1.0)/(c-1.0);
                CC = ((b - rho*sigma*i*phi1 - d) - (b - rho*sigma*i*phi1 + d)*c*Complex.Exp(-d*tau)) / (sigma*sigma) / (1.0-c*Complex.Exp(-d*tau));
            }
            else if(trap==0)
            {
                // Original Heston formulation.
                G = (1.0-g*Complex.Exp(d*tau))/(1.0-g);
                CC = ((b - rho*sigma*i*phi1 + d) - (b - rho*sigma*i*phi1 - d)*g*Complex.Exp(d*tau)) / (sigma*sigma) / (1.0-g*Complex.Exp(d*tau));
            }
            return CC;
        }
        // Elices Heston price fixed parameters
        public double ElicesPrice(string PutCall,double S,double K,double[] T,double r,double q,
            double[] param,double v0,int trap,double[] x,double[] w,double[][] paramfixed)
        {
            Complex i = new Complex(0.0,1.0);
            int NT = T.Length;
            double Mat = T[NT-1];
            
            int N = x.Length;
            double[] int1 = new double[N];
            double[] int2 = new double[N];
            Complex f1,f2;
            for(int k=0;k<=N-1;k++)
            {
                Complex phi = x[k];
                double weight = w[k];
                f2 = ElicesCF(phi  ,param,v0,T,S,r,q,trap,paramfixed);
                f1 = ElicesCF(phi-i,param,v0,T,S,r,q,trap,paramfixed)/(S*Math.Exp((r-q)*Mat));
                f2 = Complex.Exp(-i*phi*Complex.Log(K))*f2/i/phi;
                f1 = Complex.Exp(-i*phi*Complex.Log(K))*f1/i/phi;
                int2[k] = weight * f2.Real;
                int1[k] = weight * f1.Real;
            }
            // Define P1 and P2
            double PI = Math.PI;
            double P1 = 0.5 + 1.0/PI*int1.Sum();
            double P2 = 0.5 + 1.0/PI*int2.Sum();

            // The call price
            double HestonC = S*Math.Exp(-q*Mat)*P1 - K*Math.Exp(-r*Mat)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*Mat) + K*Math.Exp(-r*Mat);

            // Output the option price
            if(PutCall == "C")
                return HestonC;
            else
                return HestonP;
        }
        // Elices Heston price with NO fixed parameters
        public double ElicesPrice(string PutCall,double S,double K,double[] T,double r,double q,
            double[] param,double v0,int trap,double[] x,double[] w)
        {
            Complex i = new Complex(0.0,1.0);
            int NT = T.Length;
            double Mat = T[0];

            int N = x.Length;
            double[] int1 = new double[N];
            double[] int2 = new double[N];
            Complex f1,f2;
            for(int k=0;k<=N-1;k++)
            {
                Complex phi = x[k];
                double weight = w[k];
                f2 = ElicesCF(phi,  param,v0,T,S,r,q,trap);
                f1 = ElicesCF(phi-i,param,v0,T,S,r,q,trap)/(S*Math.Exp((r-q)*Mat));
                f2 = Complex.Exp(-i*phi*Complex.Log(K))*f2/i/phi;
                f1 = Complex.Exp(-i*phi*Complex.Log(K))*f1/i/phi;
                int2[k] = weight * f2.Real;
                int1[k] = weight * f1.Real;
            }

            // Define P1 and P2
            double PI = Math.PI;
            double P1 = 0.5 + 1.0/PI*int1.Sum();
            double P2 = 0.5 + 1.0/PI*int2.Sum();

            // The call price
            double HestonC = S*Math.Exp(-q*Mat)*P1 - K*Math.Exp(-r*Mat)*P2;

            // The put price by put-call parity
            double HestonP = HestonC - S*Math.Exp(-q*Mat) + K*Math.Exp(-r*Mat);

            // Output the option price
            if(PutCall == "C")
                return HestonC;
            else
                return HestonP;
        }
    }
}


