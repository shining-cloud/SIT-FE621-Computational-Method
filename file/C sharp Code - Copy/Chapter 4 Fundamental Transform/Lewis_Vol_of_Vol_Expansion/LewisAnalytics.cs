using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Lewis_Vol_of_Vol_Expansion
{
    class Lewis
    {
        // Ratios of Black Scholes derivatives
        public double R(double p,double q,double X,double Z,double T)
        {
            double y = 0.0;
            if(p==2 & q==0)
                y = T*(0.5*(X/Z)*(X/Z) - 0.5/Z - 1.0/8.0);
            else if(p==1 & q==1)
                y = -X/Z + 0.5;
            else if(p==1 & q==2)
                y = (X/Z)*(X/Z) - X/Z - 0.25/Z*(4.0-Z);
            else if(p==2 & q==2)
                y = T*(0.5*Math.Pow(X/Z,4) - 0.5*Math.Pow(X/Z,3) - 3.0*(X*X/Z/Z/Z) + 1.0/8.0*(X/Z/Z)*(12.0+Z) + 1.0/32.0/Z/Z*(48.0-Z*Z));
            return y;
        }
        // The "J" functions
        public double J(double rho,double theta,double k,double T,double v0,int Number)
        {
            double y = 0.0;
            double w = 0.0;
            if(Number==1)
            {
                w = (theta*T + (1.0-Math.Exp(-k*T))*(v0/k - 2.0*theta/k) - Math.Exp(-k*T)*(v0-theta)*T);
                y = w*rho/k;
            }
            else if(Number==3)
            {
                w = theta*T
	              + theta/2.0/k*(1.0-Math.Exp(-2.0*k*T)) 
            	  - 2.0*theta/k*(1.0-Math.Exp(-k*T))
	              + (v0-theta)*1.0/k*(Math.Exp(-k*T)-Math.Exp(-2.0*k*T))
	              + (v0-theta)*(-2.0*Math.Exp(-k*T)*T)
	              + (v0-theta)*1.0/k*(1.0-Math.Exp(-k*T));
                y = w/2.0/k/k;
            }
            else if(Number==4)
            {
                w = theta/k*(T*(1.0+Math.Exp(-k*T)) - 2.0/k*(1.0-Math.Exp(-k*T)))
	              + (v0-theta)/k*(1.0/k*(1.0-Math.Exp(-k*T)) - T*Math.Exp(-k*T))
	              - 0.5*T*T*Math.Exp(-k*T)*(v0-theta);
                y = w*rho*rho/k;
            }
            return y;
        }
        // Series I expansion
        public double SeriesICall(double S,double K,double rf,double q,double T,double v0,double rho,double theta,double kappa,double sigma)
        {
            // The "J" integrals
            double J1 = J(rho,theta,kappa,T,v0,1);
            double J3 = J(rho,theta,kappa,T,v0,3);
            double J4 = J(rho,theta,kappa,T,v0,4);

            // Time average of the deterministic volatility
            double v = theta + (v0-theta)*(1.0 - Math.Exp(-kappa*T))/(kappa*T);

            // X and Z required for the ratios of Black Scholes derivatives
            double X = Math.Log(S*Math.Exp((rf-q)*T)/K);
            double Z = v*T;

            // The ratios of Black Scholes derivatives
            double R20 = R(2,0,X,Z,T);
            double R11 = R(1,1,X,Z,T);
            double R12 = R(1,2,X,Z,T);
            double R22 = R(2,2,X,Z,T);

            // Black Scholes call price evaluated at v
            BlackScholes BS = new BlackScholes();
            double c = BS.BSC(S,K,rf,q,v,T);

            // Black Scholes vega evaluated at v
            double cv = BS.BSV(S,K,rf,q,v,T);

            // Series I volatility of volatility expansion for the implied volatility
            return c + sigma/T*J1*R11*cv
                  + sigma*sigma*(J3*R20/T/T + J4*R12/T + J1*J1*R22/T/T/2)*cv;
        }

        // Series II expansion
        public double[] SeriesIICall(double S,double K,double rf,double q,double T,double v0,double rho,double theta,double kappa,double sigma)
        {
            // The "J" integrals
            double J1 = J(rho,theta,kappa,T,v0,1);
            double J3 = J(rho,theta,kappa,T,v0,3);
            double J4 = J(rho,theta,kappa,T,v0,4);

            // Time average of the deterministic volatility
            double v = theta + (v0-theta)*(1.0-Math.Exp(-kappa*T))/(kappa*T);

            // X and Z required for the ratios of Black Scholes derivatives
            double X = Math.Log(S*Math.Exp((rf-q)*T)/K);
            double Z = v*T;

            // The ratios of Black Scholes derivatives
            double R20 = R(2,0,X,Z,T);
            double R11 = R(1,1,X,Z,T);
            double R12 = R(1,2,X,Z,T);
            double R22 = R(2,2,X,Z,T);

            // Series II volatility of volatility expansion for the implied volatility
            double iv = v + sigma/T*J1*R11
               + sigma*sigma*(J3*R20/T/T + J4*R12/T + 0.5/T/T*J1*J1*(R22-R11*R11*R20));

            // Series II call price
            BlackScholes BS = new BlackScholes();
            double Price = BS.BSC(S,K,rf,q,iv,T);

            // The series volatility
            double ivx = Math.Sqrt(iv);

            // Return the price and the volatility
            double[] output = new double[2] { Price,ivx };
            return output;
        }
    }
}
