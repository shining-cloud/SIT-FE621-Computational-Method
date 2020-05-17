using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Greeks_Consolidated_Integrands
{
    class HestonGreeks
    {
        public double HestonGreekConsolidated(double kappa,double theta,double lambda,double rho,double sigma,double tau,
                      double K,double S,double r,double q,double v0,int trap,double[] x,double[] w,string Greek)
        {
            Complex i = new Complex(0.0,1.0);       // Imaginary unit

            // Perform the numerical integration
            int N = x.Length;
            double[] I = new double[N];
            Complex f1=0.0,f2=0.0,J=0.0,df1=0.0,df2=0.0;
            for(int j=0;j<=N-1;j++)
            {
                double phi = x[j];
                f1 = HestonCF(phi-i,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap);
                f2 = HestonCF(phi,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap);
                if(Greek == "Price")
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (f1 - K*f2);
                else if(Greek == "Delta")
                {
                    df1 = f1 * (i*phi+1)/S;
                    df2 = f2 * i*phi/S;
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (df1 - K*df2);
                }
                else if(Greek == "Gamma")
                {
                    df1 = f1 * (-phi*phi+i*phi)/S/S;
                    df2 = f2 * (-phi*phi-i*phi)/S/S;
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (df1 - K*df2);
                }
                else if(Greek == "Theta")
                {
                    Complex[] output1 = new Complex[2];
                    Complex[] output2 = new Complex[2];
                    output1 = DiffTau(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    output2 = DiffTau(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    Complex dC1 = output1[0];
                    Complex dD1 = output1[1];
                    Complex dC2 = output2[0];
                    Complex dD2 = output2[1];
                    df1 = f1*(dC1+dD1*v0);
                    df2 = f2*(dC2+dD2*v0);
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (-r*(f1-K*f2)+df1-K*df2);
                }
                else if(Greek == "Rho")
                {
                    df1 = f1 * i*(phi-i)*tau;
                    df2 = f2 * i*phi*tau;
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (-tau*(f1-K*f2)+df1-K*df2);
                }
                else if(Greek == "Vega1")
                {
                    Complex D1 = D(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    Complex D2 = D(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    df1 = f1*D1;
                    df2 = f2*D2;
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (df1 - K*df2);
                }
                else if(Greek == "Vanna")
                {
                    Complex D1 = D(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    Complex D2 = D(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    df1 = f1 * D1*(i*phi+1)/S;
                    df2 = f2 * D2*i*phi/S;
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (df1 - K*df2);
                }
                else if(Greek == "Volga")
                {
                    Complex D1 = D(phi-i,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    Complex D2 = D(phi,kappa,theta,lambda,rho,sigma,tau,r,q,trap);
                    df1 = f1*D1;
                    df2 = f2*D2;
                    Complex d2f1 = f1*D1*D1;
                    Complex d2f2 = f2*D2*D2;
                    J = Complex.Exp(-i*phi*Math.Log(K)-r*tau)/i/phi * (4.0*v0*(d2f1-K*d2f2) + 2.0*(df1-K*df2));
                }
                I[j] = w[j] * J.Real;
            }
            double Integral = I.Sum();

            // Calculate the desired Greek
            double pi = Math.PI;
            if(Greek == "Price")
                return S*Math.Exp(-q*tau)/2.0 - K*Math.Exp(-r*tau)/2.0 + (1.0/pi)*Integral;
            else if(Greek == "Delta")
                return Math.Exp(-q*tau)/2.0 + (1.0/pi)*Integral;
            else if((Greek == "Gamma") || (Greek == "Volga"))
                return (1.0/pi)*Integral;
            else if (Greek == "Theta")
                return S/2.0*q*Math.Exp(-q*tau) - r/2.0*K*Math.Exp(-r*tau) - (1.0/pi)*Integral;
            else if (Greek == "Rho")
                return tau/2.0*K*Math.Exp(-r*tau) + (1.0/pi)*Integral;
            else if((Greek == "Vega1") || (Greek == "Vanna"))
                return (1.0/pi)*Integral*2.0*Math.Sqrt(v0);
            else
                return 0.0;
        }
        // Heston characteristic function (f2) for the Greeks
        public Complex HestonCF(Complex phi,double kappa,double theta,double lambda,double rho,double sigma,double T,
            double K,double S,double r,double q,double v0,int trap)
        {
            Complex i  = new Complex(0.0,1.0);                   // Imaginary unit
            double x = Math.Log(S);
            double a = kappa*theta;
            Complex b,u,d,g,c,D,G,C = new Complex();

            u = -0.5;
            b = kappa + lambda;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2.0) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(trap==1)
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
            return Complex.Exp(C + D*v0 + i*phi*x);
        }

        // Derivative of C and D with respect to tau
        public Complex[] DiffTau(Complex phi,double kappa,double theta,double lambda,double rho,double sigma,double tau,double r,double q,int trap)
        {
            Complex i  = new Complex(0.0,1.0);
            double a = kappa*theta;
            double u = -0.5;
            double b = kappa + lambda;
            Complex d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2.0) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            Complex g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            Complex c,D,G,C;
            if(trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                D = (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*tau))/(1.0-c*Complex.Exp(-d*tau)));
                G = (1.0 - c*Complex.Exp(-d*tau))/(1.0-c);
                C = (r-q)*i*phi*tau + a/sigma/sigma*((b - rho*sigma*i*phi - d)*tau - 2.0*Complex.Log(G));
            }
            else if(trap==0)
            {
                // Original Heston formulation.
                G = (1.0 - g*Complex.Exp(d*tau))/(1.0-g);
                C = (r-q)*i*phi*tau + a/sigma/sigma*((b - rho*sigma*i*phi + d)*tau - 2.0*Complex.Log(G));
                D = (b - rho*sigma*i*phi + d)/sigma/sigma*((1-Complex.Exp(d*tau))/(1.0-g*Complex.Exp(d*tau)));
            }

            // Derivatives of C and D w.r.t. maturity tau
            Complex dD = d*Complex.Exp(d*tau)*(b-rho*sigma*phi*i+d)*(g-1.0)/sigma/sigma/Complex.Pow(1.0-g*Complex.Exp(d*tau),2.0);
            Complex dC = (r-q)*phi*i + kappa*theta/sigma/sigma
                * ((b-rho*sigma*phi*i+d) + 2.0*g*d*Complex.Exp(d*tau)/(1.0-g*Complex.Exp(d*tau)));
            Complex[] output = new Complex[2];
            output[0] = dC;
            output[1] = dD;
            return output;
        }
        public Complex D(Complex phi,double kappa,double theta,double lambda,double rho,double sigma,double tau,double r,double q,int trap)
        {
            Complex i  = new Complex(0.0,1.0);                   // Imaginary unit
            double a = kappa*theta;
            double u = -0.5;
            double b = kappa + lambda;
            Complex d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2.0) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            Complex g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(trap==1)
            {
                Complex c = 1.0/g;
                return (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*tau))/(1.0-c*Complex.Exp(-d*tau)));
            }
            else if(trap==0)
                return (b - rho*sigma*i*phi + d)/sigma/sigma*((1.0-Complex.Exp(d*tau))/(1.0-g*Complex.Exp(d*tau)));
            else
                return 0.0;
        }
    }
}








