using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Chiarella_Ziogas_American_Call
{
    class CharFun
    {
        public Complex CZCharFun(double S0,double tau,double t,HParam param,double K,double rf,double q,Complex phi,Complex psi,int FunNum)
        {
            Complex i = new Complex(0.0,1.0);
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            double x = Math.Log(S0);

            // Parameters "a" and "b"
            double a = kappa*theta;
            double b = kappa + lambda;

            // "d" and "g" functions
            Complex d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2.0) + sigma*sigma*phi*(phi+i));
            Complex g = (b - rho*sigma*i*phi - sigma*sigma*i*psi + d) 
                      / (b - rho*sigma*i*phi - sigma*sigma*i*psi - d);

            // The components of the affine characteristic function.
            Complex G = (1.0-g*Complex.Exp(d*(tau-t)))/(1.0-g);
            Complex C = (rf-q)*i*phi*(tau-t) + a/sigma/sigma*((b - rho*sigma*i*phi + d)*(tau-t) - 2.0*Complex.Log(G));
            Complex F = (1.0-Complex.Exp(d*(tau-t)))/(1.0-g*Complex.Exp(d*(tau-t)));
            Complex D = i*psi + (b - rho*sigma*i*phi - sigma*sigma*i*psi + d)/sigma/sigma * F;

            // The characteristic function.
            Complex f2 = Complex.Exp(C + D*v0 + i*phi*x);
            if(FunNum == 2)
                return f2;
            else if(FunNum == 1)
            {
                d = Complex.Sqrt(Complex.Pow(rho*sigma*i*(phi-i) - b,2.0) + sigma*sigma*(phi-i)*phi);
                g = (b - rho*sigma*i*(phi-i) - sigma*sigma*i*psi + d) 
                  / (b - rho*sigma*i*(phi-i) - sigma*sigma*i*psi - d);

                // The components of the affine characteristic function.
                G = (1.0-g*Complex.Exp(d*(tau-t)))/(1.0-g);
                C = (rf-q)*i*(phi-i)*(tau-t) + a/sigma/sigma*((b - rho*sigma*i*(phi-i) + d)*(tau-t) - 2.0*Complex.Log(G));
                F = (1.0-Complex.Exp(d*(tau-t)))/(1.0-g*Complex.Exp(d*(tau-t)));
                D = i*psi + (b - rho*sigma*i*(phi-i) - sigma*sigma*i*psi + d)/sigma/sigma * F;

                Complex F2 = Complex.Exp(C + D*v0 + i*(phi-i)*x);
                return 1.0/S0 * Complex.Exp(-(rf-q)*(tau-t)) * F2;
            }
            else
                return 0.0;
        }
    }
}
