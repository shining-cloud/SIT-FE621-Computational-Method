using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Bounds_on_Alpha
{
    class LordKahl
    {
        // The Lord and Kahl psi function for the Heston model
        public double HestonPsi(double v,double alpha,double kappa,double theta,double lambda,double rho,double sigma,double tau,double K,double S,double r,double v0,int trap)
        {
            HestonCharFun HCF = new HestonCharFun();
            Complex i = new Complex(0.0,1.0);
            double q = 0.0;
            Complex CF = HCF.HestonCF(v-(alpha+1.0)*i,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap);
            Complex y = Complex.Exp(-r*tau) * Complex.Exp(-i*v*Complex.Log(K)) * CF / (alpha*alpha + alpha - v*v + i*v*(2.0*alpha+1.0));
            return y.Real;
        }
        // Roger Lee's "D" function
        public Complex RogerLeeD(Complex phi,double kappa,double rho,double sigma)
        {
            Complex i = new Complex(0.0,1.0);
            Complex A = (rho*sigma*phi*i - kappa);
            Complex B = sigma*sigma*(phi*i + phi*phi);
            return Complex.Sqrt(A*A + B);
        }
        // Roger Lee's "G" function
        public Complex RogerLeeG(Complex phi,double kappa,double rho,double sigma)
        {
            Complex i = new Complex(0.0,1.0);
            Complex C = (rho*sigma*phi*i - kappa);
            Complex B = sigma*sigma*(phi*i + phi*phi);
            Complex d = Complex.Sqrt(C*C + B);
            return (kappa - rho*sigma*phi*i + d) / (kappa - rho*sigma*phi*i - d);
        }
        // Roger Lee's "G*exp(d*tau)" function
        public double RogerLeeGExpD(double a,double kappa,double rho,double sigma,double tau)
        {
            Complex i = new Complex(0.0,1.0);
            Complex phi = -i*a;
            Complex C = (rho*sigma*phi*i - kappa);
            Complex B = sigma*sigma*(phi*i + phi*phi);
            Complex d = Complex.Sqrt(C*C + B);
            Complex g = (kappa - rho*sigma*phi*i + d) / (kappa - rho*sigma*phi*i - d);
            Complex E = g*Complex.Exp(d*tau);
            return (E.Real - 1.0)*(E.Real - 1.0);
        }

        // Lord and Kahl objective function to find optimal alpha
        public double LordKahlFindAlpha(double alpha,double kappa,double theta,double lambda,double rho,double sigma,
                                        double tau,double K,double S,double r,double v0)
        {
            int trap = 1;
            double PSI = HestonPsi(0.0,alpha,kappa,theta,lambda,rho,sigma,tau,K,S,r,v0,trap);
            return -alpha*Math.Log(K)  + 0.5*Math.Log(PSI*PSI);
        }
    }
}


