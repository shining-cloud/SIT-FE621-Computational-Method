using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Mikhailov_and_Nogel
{
    class CDcoefficients
    {
        // Time dependent C function
        public Complex Ct(double phi,HParam param,double r, double q,double T,OpSet settings,int Pnum,Complex C0, Complex D0)
        {
            Complex i = new Complex(0.0,1.0);                   // Imaginary unit
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double a = kappa * theta;
            int Trap = settings.trap;
            Complex b,u,d,g,G,C = new Complex();

            // Parameters "u" and "b" are different for P1 and P2
            if(Pnum==1)
            {
                u = 0.5;
                b = kappa - rho*sigma;
            }
            else
            {
                u = -0.5;
                b = kappa;
            }
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d - D0*sigma*sigma) / (b - rho*sigma*i*phi - d - D0*sigma*sigma);
            G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
            C = (r-q)*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi + d)*T - 2.0*Complex.Log(G)) + C0;

            return C;
        }
        // Time dependent D function
        public Complex Dt(double phi,HParam param,double r,double q,double T,OpSet settings,int Pnum,Complex C0,Complex D0)
        {
            Complex i = new Complex(0.0,1.0);                   // Imaginary unit
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            int Trap = settings.trap;
            Complex b,u,d,g,G,D = new Complex();

            // Parameters "u" and "b" are different for P1 and P2
            if(Pnum==1)
            {
                u = 0.5;
                b = kappa - rho*sigma;
            }
            else
            {
                u = -0.5;
                b = kappa;
            }
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d - D0*sigma*sigma) / (b - rho*sigma*i*phi - d - D0*sigma*sigma);
            G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
            D = ((b - rho*sigma*i*phi)*(1-g*Complex.Exp(d*T)) + d*(1+g*Complex.Exp(d*T))) / (sigma*sigma*(1-g*Complex.Exp(d*T)));

            return D;
        }
    }
}