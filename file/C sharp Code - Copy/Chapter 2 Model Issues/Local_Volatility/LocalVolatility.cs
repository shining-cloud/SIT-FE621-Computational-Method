using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Numerics;

namespace Local_Volatility
{
    class LocalVol
    {
        // Heston Local Volatility using Dupire's formula.
        // Uses finite differences for the required derivatives
        public double HestonLVFD(double S,double K,double T,double rf,double q,
                                 double kappa,double theta,double sigma,double v0,double rho,double lambda,
                                 int trap,double[] x,double[] w,double dt,double dK)
        {
            HestonPrice HP = new HestonPrice();
            // dC/dT by central finite difference
            double CT_1 = HP.HestonPriceGaussLaguerre("C",S,K,rf,q,T-dt,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double CT1  = HP.HestonPriceGaussLaguerre("C",S,K,rf,q,T+dt,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double dCdT =  (CT1 - CT_1) / (2.0*dt);

            // dC2/dK2 by central finite differences
            double CK_1 = HP.HestonPriceGaussLaguerre("C",S,K-dK,rf,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double CK0  = HP.HestonPriceGaussLaguerre("C",S,K,rf,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double CK1  = HP.HestonPriceGaussLaguerre("C",S,K+dK,rf,q,T,kappa,theta,sigma,v0,lambda,rho,x,w,trap);
            double dC2dK2 = (CK_1 - 2.0*CK0 + CK1) / (dK*dK);

            // Local variance
            double LocalVar = 2.0*dCdT / (K*K*dC2dK2);

            // Local volatility
            double LocalVol = Math.Sqrt(LocalVar);
            return LocalVol;
        }

        // Heston Local Volatility using Gatheral's approximation
        // Works for RHO APPROXIMATELY 1 or -1
        public double HestonLVApprox(double S,double K,double T,double kappa,double theta,double sigma,double v0,double rho)
        {
            // Modified parameters kappa' and theta'
            double kappa_ = kappa - rho*sigma/2.0;
            double theta_ = theta*kappa/kappa_;

            // wT and vT'
            double xT = Math.Log(K/S);
            double wT = (v0-theta)*(1.0-Math.Exp(-kappa*T))/kappa + theta*T;
            double vT = (v0-theta_)*Math.Exp(-kappa_*T) + theta_;

            // Integral
            double F1 = (v0-theta)/(kappa_-kappa);
            double E1 = Math.Exp((kappa_-kappa)*T) - 1.0;
            double F2 = theta/kappa_;
            double E2 = Math.Exp(kappa_*T) - 1.0;
            double Integral = Math.Exp(-kappa_*T)*(F1*E1 + F2*E2);

            // Local Variance
            double uT = vT + rho*sigma*xT/wT*Integral;

            // Local volatility
            double u = Math.Sqrt(uT);
            return u;
            if(rho == -1.0)
            {
                // Gatheral formula (4.1) when rho = -1
                double uT2 = vT - sigma*Math.Log(K/S)*(1.0-Math.Exp(-kappa_*T))/(kappa_*T);
                double y2 = Math.Sqrt(uT2);
                return y2;
            }
        }

        // Returns the integrand for the partial derivatives of P1 and P2 
        // with respect to maturity, namely dP1/dT and dP2/dT.
        public double dPjdT(int Pnum,double phi,double kappa,double theta,double lambda,double rho,double sigma,double T,double K,double S,double r,double q,double v0, int trap) 
        {
            Complex i  = new Complex(0.0,1.0);       // Imaginary unit
            Complex b,u,d,g,c,D,G,C,f,edT,dDdT,dCdT = new Complex();

            double x = Math.Log(S);
            double a = kappa*theta;

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
            if(trap==1)
            {
                // "Little Heston Trap" formulation
                c = 1.0/g;
                edT = Complex.Exp(-d*T);
                dDdT = (b - rho*sigma*i*phi - d)/sigma/sigma * (d*edT*(1.0-c*edT) - (1.0-edT)*c*d*edT)
                     / Complex.Pow(1 - c*edT,2);
                dCdT = (r-q)*i*phi + kappa*theta/sigma/sigma
                     * ((b - rho*sigma*i*phi - d) + 2.0*c*d*edT/(1.0 - c*edT));
                D = (b - rho*sigma*i*phi - d)/sigma/sigma*((1.0-Complex.Exp(-d*T))/(1.0-c*Complex.Exp(-d*T)));
                G = (1.0 - c*Complex.Exp(-d*T))/(1-c);
                C = (r-q)*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi - d)*T - 2.0*Complex.Log(G));
            }
            else
            {
                // Original Heston formulation.
                edT = Complex.Exp(d*T);
                dDdT = (b - rho*sigma*i*phi + d)/sigma/sigma * (d*edT*(g*edT-1.0) + (1.0-edT)*g*d*edT)
		             / Complex.Pow(1 - g*edT,2);
                dCdT = (r-q)*i*phi + kappa*theta/sigma/sigma
		             * ((b - rho*sigma*i*phi + d) + 2*g*d*edT/(1 - g*edT));
                G = (1.0 - g*Complex.Exp(d*T))/(1.0-g);
                C = (r-q)*i*phi*T + a/sigma/sigma*((b - rho*sigma*i*phi + d)*T - 2.0*Complex.Log(G));
                D = (b - rho*sigma*i*phi + d)/sigma/sigma*((1.0-Complex.Exp(d*T))/(1.0-g*Complex.Exp(d*T)));
            }

            // The derivative.
            Complex dfdT = Complex.Exp(C + D*v0 + i*phi*x) * (dCdT + dDdT*v0);

            // Return the real part
            Complex y = Complex.Pow(K,(-i*phi)) * dfdT / i / phi;
            return y.Real;
        }

        // Returns the integrand for the second-order partial derivative of P1 
        // with respect to strike, namely d^2P1/dK^2
        public double d2P1dK2(double phi,double Kappa,double Theta,double Lambda,double Rho,double Sigma,double Tau,double Strike,double Spot,double Rate,double Div,double V,int trap)
        {
            Complex S     = new Complex(Spot,0.0);	    // Spot Price
            Complex K     = new Complex(Strike,0.0);       // Strike Price
            Complex T     = new Complex(Tau,0.0);       // Maturity in years
            Complex r     = new Complex(Rate,0.0);       // Interest rate
            Complex q     = new Complex(Div,0.0);       // Dividend yield
            Complex i     = new Complex(0.0,1.0);       // Imaginary unit
            Complex rho   = new Complex(Rho,0.0);       // Heston parameter: correlation
            Complex kappa = new Complex(Kappa,0.0);       // Heston parameter: mean reversion speed
            Complex theta = new Complex(Theta,0.0);       // Heston parameter: mean reversion speed
            Complex lambda = new Complex(Lambda,0.0);       // Heston parameter: price of volatility risk
            Complex sigma = new Complex(Sigma,0.0);       // Heston parameter: volatility of variance
            Complex v0    = new Complex(V,0.0);        // Heston parameter: initial variance
            Complex x     = Complex.Log(S);
            Complex a     = kappa * theta;
            Complex b,u,d,g,c,D,G,C,f1,integrand = new Complex();

            // Parameters "u" and "b" for P1 
            u = 0.5;
            b = kappa + lambda - rho*sigma;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(trap==1)
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
            f1 = Complex.Exp(C + D*v0 + i*phi*x);

            // The integrand.
            integrand = (i*phi+1) * Complex.Pow(K,(-i*phi-2)) * f1;

            // The real part of the integrand
            return integrand.Real;
        }

        // Returns the integrand for the second-order partial derivative of P1 
        // with respect to strike, namely d^2P1/dK^2
        public double dP2dK2_2(double phi,double Kappa,double Theta,double Lambda,double Rho,double Sigma,double Tau,double Strike,double Spot,double Rate,double Div,double V,int trap)
        {
            Complex S     = new Complex(Spot,0.0);	    // Spot Price
            Complex K     = new Complex(Strike,0.0);       // Strike Price
            Complex T     = new Complex(Tau,0.0);       // Maturity in years
            Complex r     = new Complex(Rate,0.0);       // Interest rate
            Complex q     = new Complex(Div,0.0);       // Dividend yield
            Complex i     = new Complex(0.0,1.0);       // Imaginary unit
            Complex rho   = new Complex(Rho,0.0);       // Heston parameter: correlation
            Complex kappa = new Complex(Kappa,0.0);       // Heston parameter: mean reversion speed
            Complex theta = new Complex(Theta,0.0);       // Heston parameter: mean reversion speed
            Complex lambda = new Complex(Lambda,0.0);       // Heston parameter: price of volatility risk
            Complex sigma = new Complex(Sigma,0.0);       // Heston parameter: volatility of variance
            Complex v0    = new Complex(V,0.0);        // Heston parameter: initial variance
            Complex x     = Complex.Log(S);
            Complex a     = kappa * theta;
            Complex b,u,d,g,c,D,G,C,f2,integrand = new Complex();

            // Parameters "u" and "b" for P2 
            u = -0.5;
            b = kappa + lambda;
            d = Complex.Sqrt(Complex.Pow(rho*sigma*i*phi - b,2) - sigma*sigma*(2.0*u*i*phi - phi*phi));
            g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
            if(trap==1)
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
            f2 = Complex.Exp(C + D*v0 + i*phi*x);

            // The integrand.
            integrand = (i*phi-1) * Complex.Pow(K,(-i*phi-1)) * f2;

            // The real part of the integrand
            return integrand.Real;
        }

        // Heston Local Volatility using Dupire's equation
        public double HestonLVAnalytic(double S,double K,double T,double rf,double q,double kappa,double theta,double sigma,
                                       double lambda,double v0,double rho,double[] x,double[] w,int trap)
        {
            HestonPrice HP = new HestonPrice();
            double[] int1 = new Double[32];
            double[] int2 = new Double[32];
            double[] int3 = new Double[32];
            double[] int4 = new Double[32];
            double[] int5 = new Double[32];
            for(int k=0;k<=31;k++)
            {
                int1[k] = w[k] * dPjdT(1,x[k],kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,trap);
                int2[k] = w[k] * dPjdT(2,x[k],kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,trap);
                int3[k] = w[k] * HP.HestonProb(x[k],kappa,theta,lambda,rho,sigma,T,K,S,rf,v0,q,2,trap);
                int4[k] = w[k] * d2P1dK2(x[k],kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,trap);
                int5[k] = w[k] * dP2dK2_2(x[k],kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,trap);
            }
            double pi = Math.PI;
            double dP1dT = (1.0/pi)*int1.Sum();
            double dP2dT = (1.0/pi)*int2.Sum();
            double P2 = 0.5 + (1.0/pi)*int3.Sum();

            // dC/dT : derivative of call with respect to T
            double dCdT = S*dP1dT - K*Math.Exp(-rf*T)*(-rf*P2 + dP2dT);
            double dP1dK2 = (1.0/pi)*int4.Sum();
            double TwodP2dK2 = (1.0/pi)*int5.Sum();

            // d2C/dK2 : second derivative of call with respect to K^2
            double d2CdK2 = S*dP1dK2 - Math.Exp(-rf*T)*TwodP2dK2;

            // Local Variance
            double LocalVar = 2.0*dCdT / (K*K*d2CdK2);

            // Local Volatility
            double LocalVol = Math.Sqrt(LocalVar);
            return LocalVol;
        }
    }
}

