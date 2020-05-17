using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_LSM_Greeks
{
    class HestonGreeksClosed
    {
        // Integrands for Heston Greeks
        public double HestonGreeksProb(double phi,HParam param,OpSet settings,int Pnum,string Greek)
        {
            Complex i     = new Complex(0.0,1.0);                   // Imaginary unit
            Complex S     = new Complex(settings.S,0.0);	                // Spot Price
            Complex K     = new Complex(settings.K,0.0);               // Strike Price
            Complex T     = new Complex(settings.T,0.0);               // Maturity in years
            Complex r     = new Complex(settings.r,0.0);               // Interest rate
            Complex q     = new Complex(settings.q,0.0);               // Dividend yield
            Complex rho   = new Complex(param.rho,0.0);         // Heston parameter: correlation
            Complex kappa = new Complex(param.kappa,0.0);         // Heston parameter: mean reversion speed
            Complex theta = new Complex(param.theta,0.0);         // Heston parameter: mean reversion speed
            Complex lambda = new Complex(param.lambda,0.0);         // Heston parameter: price of volatility risk
            Complex sigma = new Complex(param.sigma,0.0);         // Heston parameter: volatility of variance
            Complex v0    = new Complex(param.v0,0.0);         // Heston parameter: initial variance
            Complex x     = Complex.Log(S);
            Complex a     = kappa * theta;
            int Trap      = settings.trap;
            Complex b,u,d,g,c,D,G,C,f = new Complex();

            // Parameters "u" and "b" are different for P1 and P2
            if(Pnum==1)
            {
                u = 0.5;
                b = kappa + lambda - rho*sigma;
            }
            else
            {
                u = -0.5;
                b = kappa + lambda;
            }
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(Trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                D = (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*T))/(1.0-c*Complex.Exp(-d*T)));
                G = (1.0 - c*Complex.Exp(-d*T))/(1-c);
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
            f = Complex.Exp(C + D*v0 + i*phi*x);

            // Calculate the integrands for the various Greeks
            Complex y = new Complex(0.0,0.0);
            Complex dD = new Complex(0.0,0.0);
            Complex dC = new Complex(0.0,0.0);
            Complex df = new Complex(0.0,0.0);

            if((Greek == "Delta") | (Greek == "Rho"))
                y = Complex.Exp(-i*phi*Complex.Log(K))*f/i/phi;
            else if(Greek == "Gamma")
                // der(P1^2)/der(S^2)
                y = Complex.Exp(-i*phi*Complex.Log(K))*f;
            else if(Greek == "Theta")
            {
                // dP/dtau
                dD = d*Complex.Exp(d*T)*(b-rho*sigma*phi*i+d)*(g-1.0)/sigma/sigma/Complex.Pow(1.0-g*Complex.Exp(d*T),2);
                dC = (r-q)*phi*i + kappa*theta/sigma/sigma * ((b-rho*sigma*phi*i+d) + 2.0*g*d*Complex.Exp(d*T)/(1.0-g*Complex.Exp(d*T)));
                df = f*(dC + dD*v0);
                y = Complex.Exp(-i*phi*Complex.Log(K))*df/i/phi;
            }
            else if(Greek == "Vega1")
                // dP/dv0
                y = Complex.Exp(-i*phi*Complex.Log(K))*f*D/i/phi;
            else if(Greek == "Vega2")
            {
                if(Trap == 1)
                    dC = (r-q)*i*phi*T + kappa/sigma/sigma*((b - rho*sigma*i*phi - d)*T - 2.0*Complex.Log(G));
                else if(Trap == 0)
                    dC = (r-q)*i*phi*T + kappa/sigma/sigma*((b - rho*sigma*i*phi + d)*T - 2.0*Complex.Log(G));
                df =f*dC;
                y = Complex.Exp(-i*phi*Complex.Log(K))*df/i/phi;
            }
            else if(Greek =="Volga")
                // dP2/dv02
                y = Complex.Exp(-i*phi*Complex.Log(K))*f*D*D/i/phi;
            return y.Real;
        }

        // Heston Closed Form Greeks
        public double HestonGreeks(HParam param,OpSet settings,double[] X,double[] W,string Greek)
        {
            int N = X.Length;
            double P1,P2,dC2,dP1,dP2;
            double pi = Math.PI;
            double r = settings.r;
            double q = settings.q;
            double T = settings.T;
            double K = settings.K;
            double S = settings.S;
            double theta = param.theta;
            string PutCall = settings.PutCall;
            double v0 = param.v0;
            double[] int1 = new double[N];
            double[] int2 = new double[N];
            double[] der1 = new double[N];
            double[] der2 = new double[N];
            double y = 0.0;

            if(Greek == "Delta")
            {
                for(int k=0;k<=N-1;k++)
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Delta");
                P1 = 0.5 + 1.0/pi*int1.Sum();
                if(PutCall == "C")
                    y = Math.Exp(-q*T)*P1;
                else
                    y = Math.Exp(-q*T)*(P1 - 1.0);
            }
            else if(Greek == "Gamma")
            {
                for(int k=0;k<=N-1;k++)
                    // der(P1^2)/der(S^2)
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Gamma");
                y = Math.Exp(-q*T)/pi/S*int1.Sum();
            }
            else if(Greek == "Rho")
            {
                for(int k=0;k<=N-1;k++)
                    int2[k] = W[k] * HestonGreeksProb(X[k],param,settings,2,"Rho");
                P2 = 0.5 + 1.0/pi*int2.Sum();
                if(PutCall == "C")
                    y = K*Math.Exp(-r*T)*T*P2;
                else
                    y = K*Math.Exp(-r*T)*T*(P2-1.0);
            }
            else if(Greek == "Theta")
            {
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Delta");
                    int2[k] = W[k] * HestonGreeksProb(X[k],param,settings,2,"Delta");
                    der1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Theta");
                    der2[k] = W[k] * HestonGreeksProb(X[k],param,settings,2,"Theta");
                }
                P1 = 0.5 + 1.0/pi*int1.Sum();
                P2 = 0.5 + 1.0/pi*int2.Sum();
                dP1 = 1.0/pi*der1.Sum();
                dP2 = 1.0/pi*der2.Sum();
                double theta2 = -S*Math.Exp(-q*T)*(-q*P1 + dP1) + K*Math.Exp(-r*T)*(dP2 - r*P2);
                if(PutCall == "C")
                    y = theta2;
                else if(PutCall == "P")
                    y = theta2 + K*r*Math.Exp(-r*T) - q*S*Math.Exp(-q*T);
            }
            else if(Greek == "Vega1")
            {
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Vega1");
                    int2[k] = W[k] * HestonGreeksProb(X[k],param,settings,2,"Vega1");
                }
                dP1 = 1.0/pi*int1.Sum()*2.0*Math.Sqrt(v0);
                dP2 = 1.0/pi*int2.Sum()*2.0*Math.Sqrt(v0);
                y = S*Math.Exp(-q*T)*dP1 - K*Math.Exp(-r*T)*dP2;
            }
            else if(Greek == "Vega2")
            {
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Vega2");
                    int2[k] = W[k] * HestonGreeksProb(X[k],param,settings,2,"Vega2");
                }
                dP1 = 1.0/pi*int1.Sum()*2.0*Math.Sqrt(theta);
                dP2 = 1.0/pi*int2.Sum()*2.0*Math.Sqrt(theta);
                y = S*Math.Exp(-q*T)*dP1 - K*Math.Exp(-r*T)*dP2;
            }
            else if(Greek == "Vanna")
            {
                for(int k=0;k<=N-1;k++)
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Vega1");
                dP1 = 1.0/pi*int1.Sum()*2.0*Math.Sqrt(v0);
                y = Math.Exp(-q*T)*dP1;
            }
            else if(Greek == "Volga")
            {
                // Calculate Vega1
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Vega1");
                    int2[k] = W[k] * HestonGreeksProb(X[k],param,settings,2,"Vega1");
                }
                dP1 = 1.0/pi*int1.Sum()*2.0*Math.Sqrt(v0);
                dP2 = 1.0/pi*int2.Sum()*2.0*Math.Sqrt(v0);
                double Vega1 = S*Math.Exp(-q*T)*dP1 - K*Math.Exp(-r*T)*dP2;

                // Calculate second-order derivative of the call price w.r.t. v0
                for(int k=0;k<=N-1;k++)
                {
                    int1[k] = W[k] * HestonGreeksProb(X[k],param,settings,1,"Volga");
                    int2[k] = W[k] * HestonGreeksProb(X[k],param,settings,2,"Volga");
                }
                dP1 = 1.0/pi*int1.Sum();
                dP2 = 1.0/pi*int2.Sum();
                dC2 = S*Math.Exp(-q*T)*dP1 - K*Math.Exp(-r*T)*dP2;

                // Calculate Volga
                y = 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1/4.0/v0);
            }
        return y;
        }
    }
}

