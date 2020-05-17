using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Chiarella_Ziogas_American_Call
{
    class findBfunctions
    {
        public BVec findB(double tau,HParam param,double K,double r,double q,double V00,double V10,double b00,double b10,
            double[] xs,double[] ws,double[] xt,double[] wt,int Nt,int Ntau,double tol0,double tol1,double Ntol,
            double a,double b,double c,double d,string DoubleType)
        {
            double[] v0 = new double[Ntau];
            double[] v1 = new double[Ntau];
            double[] b0 = new double[Ntau];
            double[] b1 = new double[Ntau];

            // Initial values
            v0[0] = V00;
            v1[0] = V10;
            b0[0] = b00;
            b1[0] = b10;

            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0t   = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            double dtau = tau/Ntau;
            Console.WriteLine("  n     mat     v0(n)    v1(n)    b0(n)    b1(n)   NumSteps");
            Console.WriteLine("------------------------------------------------------------");
            Console.WriteLine("{0,3:F0} {1,8:F4} {2,8:F4} {3,8:F4} {4,8:F4} {5,8:F4}",1,dtau,v0[0],v1[0],b0[0],b1[0]);

            for(int n=1;n<=Ntau-1;n++)
            {
                // Maturity incremenent
                double mat = Convert.ToDouble(n)*dtau;
                // Variances
                double Evt = EV(v0t,theta,kappa,tau,mat);
                v0[n] = Evt + sigma/Math.Abs(kappa)*Math.Sqrt(kappa*theta/2.0);
                v1[n] = Evt - sigma/Math.Abs(kappa)*Math.Sqrt(kappa*theta/2.0);
                // Parameters with v1
                HParam params1;
                params1.kappa = kappa;
                params1.theta = theta;
                params1.sigma = sigma;
                params1.v0 = v1[n];
                params1.rho = rho;
                params1.lambda = lambda;
                // Parameters with v0
                HParam params0;
                params0.kappa = kappa;
                params0.theta = theta;
                params0.sigma = sigma;
                params0.v0 = v0[n];
                params0.rho = rho;
                params0.lambda = lambda;
                // Starting values for Newton's method
                double b0k_ = b0[n-1];
                double b1k_ = b1[n-1];
                // Loop through until Newton's method returns bk and bk(-1) to tolerance
                int counter = 0;
                double diff0 = 1.1*tol0;
                double diff1 = 1.1*tol1;
                NewtonMethod NM = new NewtonMethod();
                while((diff0>tol0) && (diff1>tol1))
                {
                    counter += 1;
                    // Newton's method
                    double b1k = NM.CZNewton(b1k_,v0[n],v1[n],mat,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k_,1,Ntol,a,b,c,d,DoubleType);
                    double b0k = NM.CZNewton(b0k_,v0[n],v1[n],mat,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k ,0,Ntol,a,b,d,c,DoubleType);
                    //double b1k = NM.CZNewton(b1k_,v0[n],v1[n],tau,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k_,1,Ntol,a,b,c,d,DoubleType);
                    //double b0k = NM.CZNewton(b0k_,v0[n],v1[n],tau,kappa,theta,lambda,rho,sigma,K,r,q,xs,ws,xt,wt,Nt,b0k_,b1k ,0,Ntol,a,b,d,c,DoubleType);
                    b0[n] = b0k;
                    b1[n] = b1k;
                    diff0 = Math.Abs(b0k_ - b0k);
                    diff1 = Math.Abs(b1k_ - b1k);
                    // Update the values
                    b0k_ = b0k;
                    b1k_ = b1k;
                }
                Console.WriteLine("{0,3:F0} {1,8:F4} {2,8:F4} {3,8:F4} {4,8:F4} {5,8:F4} {6,5:F0}",n,mat,v0[n],v1[n],b0[n],b1[n],counter);
            }
            BVec bvector;
            bvector.B0 = b0;
            bvector.B1 = b1;
            return bvector;
        }
        // Function for CIR expected value E[V(t)|V(s)] for t>s
        public double EV(double Vs,double theta,double kappa,double t,double s)
        {
            return theta + (Vs - theta)*Math.Exp(-kappa*(t-s));
        }
    }
}
