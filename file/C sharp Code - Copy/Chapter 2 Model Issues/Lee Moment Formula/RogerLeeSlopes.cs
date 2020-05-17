using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lee_Moment_Formula
{
    class LeeMoment
    {
        public double[] FindLeeBounds(double S,double r,double div,double tau,HParam param,int trap,int LoLimit,int HiLimit)
        {
            // Find the upper moment bound using Andersen and Piterbarg (2007) formula
            double[] W = new double[HiLimit];
            for(int k=0;k<=HiLimit-1;k++)
                W[k] = Convert.ToDouble(k+1);

            // Infinity
            double zero = 0.0;
            double inf = 1.0 / zero;

            // Moment explosion setting
            double lambda = 1.0;

            int j = 0;
            double e = 0.0;
            double start = 0.0;
            double T = inf;
            double UpperAP = 0.0;
            for(int k=0;k<=15;k++)
            {
                j = 0;
                T = inf;
                while(T == inf)
                {
                    j = j+1;
                    T = MomentExplode(W[j],lambda,param.sigma,param.kappa,param.rho);
                }
                start = W[j-1];
                e = 1.0 / Math.Pow(10,Convert.ToDouble(k+1));
                W = new double[21];
                W[0] = start;
                for(int s=1;s<=20;s++)
                    W[s] = W[s-1] + e;
                UpperAP = W[j];
            }

            // Find the lower moment bound using Andersen and Piterbarg (2007) formula
            int MM = -LoLimit+1;
            W = new double[Math.Abs(LoLimit)+1];
            for(int k=0;k<=Math.Abs(LoLimit);k++)
                W[k] = Convert.ToDouble(LoLimit+k);

            double LowerAP = 0.0;
            for(int k=0;k<=15;k++)
            {
                j = 0;
                T = 0.0;
                while(T < inf)
                {
                    j = j+1;
                    T = MomentExplode(W[j],lambda,param.sigma,param.kappa,param.rho);
                }
                start = W[j-1];
                e = 1.0 / Math.Pow(10,Convert.ToDouble(k+1));
                W = new double[21];
                W[0] = start;
                for(int s=1;s<=20;s++)
                    W[s] = W[s-1] + e;
                LowerAP = W[j];
            }

            // Roger Lee moment formulas
            double p =  UpperAP - 1;
            double q = -LowerAP;
            double bR = 2.0 - 4.0*(Math.Sqrt(p*p + p) - p);
            double bL = 2.0 - 4.0*(Math.Sqrt(q*q + q) - q);

            double[] output = new double[4];
            output[0] = bR;
            output[1] = bL;
            output[2] = LowerAP;
            output[3] = UpperAP;

            return output;
        }
        // Andersen and Piterbarg function to find T* the time of moment explosion
        public double MomentExplode(double w,double lambda,double sigma,double kappa,double rho)
        {
            double g,beta,PI;
            double k = lambda*lambda*w*(w-1)*0.5;
            double b = 2.0*k/sigma/sigma;
            double a = 2.0*(rho*sigma*lambda*w - kappa)/sigma/sigma;
            double D = a*a - 4.0*b;
            float zero = 0;
            float inf = 1 / zero;
            double T = inf;
            if(D>=0.0 & a<0.0)
                T = inf;
            else if(D>=0.0 & a>0.0)
            {
                g = Math.Pow(D,0.5)/2.0;
                T = Math.Log((a/2.0+g)/(a/2.0-g))/g/sigma/sigma;
            }
            else if(D<0.0)
            {
                beta = Math.Pow(-D,0.5)/2;
                if(a<0.0)
                    PI = Math.PI;
                else
                    PI = 0;
                T = 2.0*(PI + Math.Atan(2.0*beta/a))/beta/sigma/sigma;
            }
            return T;
        }
    }
}

