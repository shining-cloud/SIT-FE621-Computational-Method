using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Carr_and_Madan_FFT_or_FRFT_Greeks
{
    class HestonGreeks
    {
        public double CarrMadanGreeks(double alpha,double tau,HParam param,OpSet Opsettings,string Greek,double[] x,double[] w)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit
            double q = 0.0;
            int N = x.Length;

            // Perform the numerical integration
            Complex J = 0.0;
            double[] I = new double[N];

            for(int j=0;j<=N-1;j++)
            {
                double u = x[j];
                J = w[j] * Complex.Exp(-i*u*Math.Log(Opsettings.K)) * Math.Exp(-Opsettings.r*Opsettings.T)
                  * HestonCFGreek(u-(alpha+1)*i,param,Opsettings,Greek)
                  / (alpha*alpha + alpha - u*u + i*(2.0*alpha+1.0)*u);
                I[j] = J.Real;
            }

            // Calculate the desired Greek
            double pi = Math.PI;
            if((Greek=="Delta") || (Greek=="Gamma") || (Greek=="Rho") ||
               (Greek=="Theta") || (Greek=="Volga") || (Greek=="Price"))
                return Math.Exp(-alpha*Math.Log(Opsettings.K))*I.Sum()/pi;
            else if((Greek=="Vega1") || (Greek=="Vanna"))
                return Math.Exp(-alpha*Math.Log(Opsettings.K))*I.Sum()/pi * 2.0*Math.Sqrt(param.v0);
            else
                return 0.0;
        }


        // Heston characteristic function (f2) for the Greeks
        public Complex HestonCFGreek(Complex phi,HParam param,OpSet settings,string Greek)
        {
            Complex i  = new Complex(0.0,1.0);                   // Imaginary unit
            double S = settings.S;
            double K = settings.K;
            double T = settings.T;
            double r = settings.r;
            double q = settings.q;
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0 = param.v0;
            double rho = param.rho;
            double lambda = param.lambda;
            double x = Math.Log(S);
            double a = kappa*theta;
            int Trap = settings.trap;
            Complex b,u,d,g,c,D,G,C = new Complex();

            u = -0.5;
            b = kappa + lambda;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2.0) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                D = (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*T))/(1.0-c*Complex.Exp(-d*T)));
                G = (1.0 - c*Complex.Exp(-d*T))/(1.0-c);
                C = (r-q)*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi - d)*T - 2.0*Complex.Log(G));
            }
            else
            {
                // Original Heston formulation.
                G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
                C = (r-q)*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi + d)*T - 2.0*Complex.Log(G));
                D = (b - rho*sigma*i*phi + d)/sigma/sigma*((1.0-Complex.Exp(d*T))/(1.0-g*Complex.Exp(d*T)));
            }

            // The characteristic function.
            Complex f = Complex.Exp(C + D*v0 + i*phi*x);

            // Return the integrand of choice
            if(Greek == "Price")
                return f;
            else if(Greek=="Delta")
                return f*i*phi/S;
            else if(Greek=="Gamma")
                return i*phi*f*(i*phi-1.0)/S/S;
            else if(Greek=="Theta")
            {
                Complex dDdT = d*Complex.Exp(d*T)*(b-rho*sigma*phi*i+d)*(g-1.0)/sigma/sigma/Complex.Pow(1.0-g*Complex.Exp(d*T),2.0);
                Complex dCdT = r*phi*i + kappa*theta/sigma/sigma * ((b-rho*sigma*phi*i+d) + 2.0*g*d*Complex.Exp(d*T)/(1.0-g*Complex.Exp(d*T)));
                return r*f - f*(dCdT + dDdT*v0);
            }
            else if(Greek=="Rho")
            {
                Complex dCdr = i*phi*T;
                return (f*dCdr - T*f);
            }
            else if(Greek=="Vega1")
                return f*D;
            else if(Greek=="Vanna")
                return f*i*phi*D/S;
            else if(Greek=="Volga")
                return 2.0*D*f*(2.0*D*v0 + 1.0);
            else
                return 0.0;
        }
    }
}


