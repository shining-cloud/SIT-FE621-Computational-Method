using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Double_Heston_Greeks
{
    class DHGreeks
    {
        public double DoubleHestonGreeks(DHParam param,OpSet settings,double[] X,double[] W,string Greek)
        {
            DHCharFun DH = new DHCharFun();

            int N = X.Length;
            double S = settings.S;
            double K = settings.K;
            double r = settings.r;
            double q = settings.q;
            double T = settings.T;
            int trap = settings.trap;
            string PutCall = settings.PutCall;
            Complex i = new Complex(0.0,1.0);
            Complex u = 0.0;
            Complex f1,f2,df1,df2;
            double[] Int1 = new double[N];
            double[] Int2 = new double[N];
            double[] dInt1 = new double[N];
            double[] dInt2 = new double[N];
            Complex I1=0.0,I2=0.0,B1=0.0,B2=0.0,dI1=0.0,dI2=0.0;
            double pi = Math.PI;
            double v0,v01,v02,P1,P2,dP1,dP2,dC;
            int v0choice;
            for(int k=0;k<=N-1;k++)
            {
                u = X[k];
                f2 = DH.DoubleHestonCF(u,param,settings);
                f1 = DH.DoubleHestonCF(u-i,param,settings);
                if(Greek=="Price")
                {
                    I2 = Complex.Exp(-i*u*Complex.Log(K))/i/u*f2;
                    I1 = Complex.Exp(-i*u*Complex.Log(K))/i/u/Math.Exp((r-q)*T)*f1/S;
                }
                else if(Greek=="Delta")
                    I1 = Complex.Exp(-i*u*Complex.Log(K))/i/u/Math.Exp((r-q)*T)*f1/S;
                else if(Greek=="Gamma")
                    I1 = Complex.Exp(-i*u*Complex.Log(K))/Math.Exp((r-q)*T)*f1/S/S;
                else if(Greek=="Rho")
                    I2 = Complex.Exp(-i*u*Complex.Log(K))/i/u*f2;
                else if((Greek == "Vega11") || (Greek =="Vega12"))
                {
                    if(Greek =="Vega11")
                        v0choice = 1;
                    else
                        v0choice = 2;
                    B1 = B(u-i,param,T,trap,v0choice);
                    B2 = B(u,param,T,trap,v0choice);
                    df1 = f1*B1;
                    df2 = f2*B2;
                    dI2 = Complex.Exp(-i*u*Complex.Log(K))/i/u*df2;
                    dI1 = Complex.Exp(-i*u*Complex.Log(K))/i/u/Math.Exp((r-q)*T)*df1/S;
                }
                else if((Greek == "Vanna1") || (Greek =="Vanna2"))
                {
                    if(Greek =="Vanna1")
                        v0choice = 1;
                    else
                        v0choice = 2;
                    B1 = B(u-i,param,T,trap,v0choice);
                    df1 = f1*B1;
                    dI1 = Complex.Exp(-i*u*Complex.Log(K))/i/u/Math.Exp((r-q)*T)*df1/S;
                }
                else if((Greek =="Volga1") || (Greek == "Volga2"))
                {
                    if(Greek =="Volga1")
                    {
                        v0choice = 1;
                        v0 = param.v01;
                    }
                    else
                    {
                        v0choice = 2;
                        v0 = param.v02;
                    }
                    B1 = B(u-i,param,T,trap,v0choice);
                    B2 = B(u,param,T,trap,v0choice);
                    df1 = 2.0*B1*f1*(2.0*B1*v0 + 1.0);
                    df2 = 2.0*B2*f2*(2.0*B2*v0 + 1.0);
                    dI1 = Complex.Exp(-i*u*Complex.Log(K))/i/u/Complex.Exp((r-q)*T)*df1/S;
                    dI2 = Complex.Exp(-i*u*Complex.Log(K))/i/u*df2;
                }
                else if(Greek =="Theta")
                {
                    I2 = Complex.Exp(-i*u*Complex.Log(K))/i/u*f2;
                    I1 = Complex.Exp(-i*u*Complex.Log(K))/i/u/Math.Exp((r-q)*T)*f1/S;
                    v01 = param.v01;
                    v02 = param.v02;
                    Complex[] output = new Complex[3];
                    output = DiffTau(u-i,param,T,r,q);
                    Complex dA1 = output[0];
                    Complex dB11 = output[1];
                    Complex dB21 = output[2];
                    df1 = f1*(dA1 + dB11*v01 + dB21*v02) - (r-q)*f1;
                    output = DiffTau(u,param,T,r,q);
                    Complex dA2 = output[0];
                    Complex dB12 = output[1];
                    Complex dB22 = output[2];
                    df2 = f2*(dA2 + dB12*v01 + dB22*v02);
                    dI2 = Complex.Exp(-i*u*Math.Log(K))/i/u*df2;
                    dI1 = Complex.Exp(-i*u*Math.Log(K))/i/u/Math.Exp((r-q)*T)*df1/S;
                }
                Int2[k] = W[k] * I2.Real;
                Int1[k] = W[k] * I1.Real;
                dInt2[k] = W[k] * dI2.Real;
                dInt1[k] = W[k] * dI1.Real;
            }
            // Return the price or the Greek
            if(Greek=="Price")
            {
                P1 = 0.5 + 1.0/pi*Int1.Sum();
                P2 = 0.5 + 1.0/pi*Int2.Sum();
                return S*Math.Exp(-q*T)*P1 - K*Math.Exp(-r*T)*P2;
            }
            else if(Greek=="Delta")
            {
                P1 = 0.5 + 1.0/pi*Int1.Sum();
                return Math.Exp(-q*T)*P1;
            }
            else if(Greek=="Gamma")
                return Math.Exp(-q*T)*Int1.Sum()/pi;
            else if (Greek=="Rho")
            {
                P2 = 0.5 + 1.0/pi*Int2.Sum();
                return K*T*Math.Exp(-r*T)*P2;
            }
            else if((Greek == "Vega11") || (Greek =="Vega12"))
            {
                if(Greek == "Vega11")
                    v0 = param.v01;
                else
                    v0 = param.v02;
                dP1 = 1.0/pi*dInt1.Sum();
                dP2 = 1.0/pi*dInt2.Sum();
                dC = S*Math.Exp(-q*T)*dP1 - K*Math.Exp(-r*T)*dP2;
                return dC*2.0*Math.Sqrt(v0);
            }
            else if((Greek == "Vanna1") || (Greek =="Vanna2"))
            {
                if(Greek == "Vanna1")
                    v0 = param.v01;
                else
                    v0 = param.v02;
                dP1 = 1.0/pi*dInt1.Sum();
                return 2.0*Math.Exp(-q*T)*Math.Sqrt(v0)*dP1;
            }
            else if((Greek =="Volga1") || (Greek == "Volga2"))
            {
                dP1 = 1.0/pi*dInt1.Sum();
                dP2 = 1.0/pi*dInt2.Sum();
                return S*Math.Exp(-q*T)*dP1 - K*Math.Exp(-r*T)*dP2;
            }
            else if(Greek=="Theta")
            {
                P1 = 0.5 + 1.0/pi*Int1.Sum();
                P2 = 0.5 + 1.0/pi*Int2.Sum();
                dP1 = 1.0/pi*dInt1.Sum();
                dP2 = 1.0/pi*dInt2.Sum();
                return -S*Math.Exp(-q*T)*(-q*P1+dP1) + K*Math.Exp(-r*T)*(-r*P2+dP2);
            }
            else
                return 0.0;
        }
        // Returns the "B" coefficient of the DH CF
        public Complex B(Complex phi,DHParam param,double tau,int trap,int j)
        {
            Complex i = new Complex(0.0,1.0);

            // First set of parameters
            double kappa1 = param.kappa1;
            double theta1 = param.theta1;
            double sigma1 = param.sigma1;
            double v01    = param.v01;
            double rho1   = param.rho1;

            // Second set of parameters
            double kappa2 = param.kappa2;
            double theta2 = param.theta2;
            double sigma2 = param.sigma2;
            double v02    = param.v02;
            double rho2   = param.rho2;

            Complex d1,d2,g1,g2,B1,B2;
            if(trap==1)
            {
                d1 = Complex.Sqrt(Complex.Pow(kappa1-rho1*sigma1*i*phi,2.0) + sigma1*sigma1*phi*(phi+i));
                d2 = Complex.Sqrt(Complex.Pow(kappa2-rho1*sigma2*i*phi,2.0) + sigma2*sigma2*phi*(phi+i));
                g1 = (kappa1-rho1*sigma1*phi*i-d1) / (kappa1-rho1*sigma1*phi*i+d1);
                g2 = (kappa2-rho2*sigma2*phi*i-d2) / (kappa2-rho2*sigma2*phi*i+d2);
                B1 = (kappa1-rho1*sigma1*phi*i-d1)*(1.0-Complex.Exp(-d1*tau)) / sigma1/sigma1 / (1.0-g1*Complex.Exp(-d1*tau));
                B2 = (kappa2-rho2*sigma2*phi*i-d2)*(1.0-Complex.Exp(-d2*tau)) / sigma2/sigma2 / (1.0-g2*Complex.Exp(-d2*tau));
            }
            else
            {
                d1 = Complex.Sqrt(Complex.Pow(kappa1-rho1*sigma1*phi*i,2.0) + sigma1*sigma1*(phi*i+phi*phi));
                d2 = Complex.Sqrt(Complex.Pow(kappa2-rho2*sigma2*phi*i,2.0) + sigma2*sigma2*(phi*i+phi*phi));
                g1 = (kappa1-rho1*sigma1*phi*i+d1)/(kappa1-rho1*sigma1*phi*i-d1);
                g2 = (kappa2-rho2*sigma2*phi*i+d2)/(kappa2-rho2*sigma2*phi*i-d2);
                B1 = (kappa1-rho1*sigma1*phi*i+d1)*(1.0-Complex.Exp(d1*tau))/sigma1/sigma1/(1.0-g1*Complex.Exp(d1*tau));
                B2 = (kappa2-rho2*sigma2*phi*i+d2)*(1.0-Complex.Exp(d2*tau))/sigma2/sigma2/(1.0-g2*Complex.Exp(d2*tau));
            }

            if(j==1)
                return B1;
            else
                return B2;
        }
        // Derivatives of A, B1, and B2 w.r.t. tau
        public Complex[] DiffTau(Complex phi,DHParam param,double tau,double r,double q)
        {
            Complex i = new Complex(0.0,1.0);

            // First set of parameters
            double kappa1 = param.kappa1;
            double theta1 = param.theta1;
            double sigma1 = param.sigma1;
            double v01    = param.v01;
            double rho1   = param.rho1;

            // Second set of parameters
            double kappa2 = param.kappa2;
            double theta2 = param.theta2;
            double sigma2 = param.sigma2;
            double v02    = param.v02;
            double rho2   = param.rho2;

            Complex d1 = Complex.Sqrt(Complex.Pow(kappa1-rho1*sigma1*phi*i,2.0) + sigma1*sigma1*(phi*i+phi*phi));
            Complex d2 = Complex.Sqrt(Complex.Pow(kappa2-rho2*sigma2*phi*i,2.0) + sigma2*sigma2*(phi*i+phi*phi));
            Complex g1 = (kappa1-rho1*sigma1*phi*i+d1)/(kappa1-rho1*sigma1*phi*i-d1);
            Complex g2 = (kappa2-rho2*sigma2*phi*i+d2)/(kappa2-rho2*sigma2*phi*i-d2);

            // The derivatives
            Complex dC1 = kappa1*theta1/sigma1/sigma1 * ((kappa1-rho1*sigma1*phi*i+d1) + 2.0*g1*d1*Complex.Exp(d1*tau)/(1.0-g1*Complex.Exp(d1*tau)));
            Complex dC2 = kappa2*theta2/sigma2/sigma2 * ((kappa2-rho2*sigma2*phi*i+d2) + 2.0*g2*d2*Complex.Exp(d2*tau)/(1.0-g2*Complex.Exp(d2*tau)));

            Complex dA = (r-q)*phi*i + dC1 + dC2;
            Complex dB1 = d1*Complex.Exp(d1*tau)*(kappa1-rho1*sigma1*phi*i+d1)*(g1-1.0)/sigma1/sigma1/Complex.Pow(1.0-g1*Complex.Exp(d1*tau),2.0);
            Complex dB2 = d2*Complex.Exp(d2*tau)*(kappa2-rho2*sigma2*phi*i+d2)*(g2-1.0)/sigma2/sigma2/Complex.Pow(1.0-g2*Complex.Exp(d2*tau),2.0);
            Complex[] output = new Complex[3];
            output[0] = dA;
            output[1] = dB1;
            output[2] = dB2;
            return output;
        }
    }
}

