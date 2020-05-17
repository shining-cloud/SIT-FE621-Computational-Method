using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Gauthier_Starting_Values
{
    class Gauthier
    {
        // Get closed-form Gauthier starting values
        public double[] GetGauthierValues(double kappa,double theta,double v0,double S,double K1,double K2,double Put1,double Put2,double tau,double rf,double q)
        {
            double[] Co1 = GauthierCoefficients(kappa,theta,v0,S,K1,tau,rf,q);
            double[] Co2 = GauthierCoefficients(kappa,theta,v0,S,K2,tau,rf,q);

            double A1 = Co1[0];
            double B1 = Co1[1];
            double C1 = Co1[2];
            double D1 = Co1[3];

            double A2 = Co2[0];
            double B2 = Co2[1];
            double C2 = Co2[2];
            double D2 = Co2[3];


            // Settings for the Nelder Mead algorithm
            int N = 2;						// Number of Heston parameters 
            int NumIters = 1;			    // First Iteration
            int MaxIters = 1000;		    // Maximum number of iterations
            double Tolerance = 1e-5;		// Tolerance on best and worst function values

            // Bounds on the estimates for sigma and rho
            double e = 1e-5;
            double[] lb = new double[] {e  , -.999};  // Lower bound on the estimates
            double[] ub = new double[] {3.0,  .999};  // Upper bound on the estimates

            // Starting values for {sigma}, {rho} for Neldger Mean
            double[,] s = new double[2,3] {{0.5,0.4,0.1}, {-0.3,-0.6,-0.75} };

            double[] Coeff1 = new double[4] {A1,B1,C1,D1};
            double[] Coeff2 = new double[4] {A2,B2,C2,D2};

            // Find the coefficients
            NelderMeadAlgo NM = new NelderMeadAlgo();
            ObjectiveFunction OF = new ObjectiveFunction();
            double[] Result = NM.NelderMead(OF.f,N,NumIters,MaxIters,Tolerance,s,Coeff1,Coeff2,Put1,Put2);
            return Result;

        }
        // Calculate the Gauthier coefficients
        public double[] GauthierCoefficients(double kappa,double theta,double v0,double S,double K,double T,double rf,double q)
        {
            BlackScholesPrice BS = new BlackScholesPrice();
            double x = Math.Log(S);
            double kT = kappa*T;

            double m0 = Math.Exp(-kT)*(Math.Exp(kT)-1)/kappa;
            double m1 = T - m0;
            double varT = m0*v0 + m1*theta;

            double p0 = Math.Exp(-kT)*(Math.Exp(kT) - kT - 1)/kappa/kappa;
            double p1 = Math.Exp(-kT)*(Math.Exp(kT)*(kT-2) + kT + 2)/kappa/kappa;

            double q0 = Math.Exp(-kT)*(2*Math.Exp(kT) - kT*(kT+2) - 2)/2/kappa/kappa/kappa;
            double q1 = Math.Exp(-kT)*(2*Math.Exp(kT)*(kT-3) + kT*(kT+4) + 6)/2/kappa/kappa/kappa;

            double r0 = Math.Exp(-2*kT)*(2*Math.Exp(2*kT) - 4*Math.Exp(kT)*kT - 2)/4/kappa/kappa/kappa;
            double r1 = Math.Exp(-2*kT)*(Math.Exp(2*kT)*(2*kT-5) + 4*Math.Exp(kT)*(kT+1) + 1)/4/kappa/kappa/kappa;

            double[] BSderiv = BS.BlackScholesDerivatives(kappa,theta,v0,S,K,T,rf,q);
            double P11 = BSderiv[0];
            double P21 = BSderiv[1];
            double P02 = BSderiv[2];
            double P22 = BSderiv[3];
            double y = varT;

            // Black Scholes Put Price
            double g = Math.Pow(y,-0.5)*(-x + Math.Log(K) - (rf-q)*T) - 0.5*Math.Sqrt(y);
            double f = Math.Pow(y,-0.5)*(-x + Math.Log(K) - (rf-q)*T) + 0.5*Math.Sqrt(y);
            double BSPut = K*Math.Exp(-rf*T)*BS.NormCDF(f) - S*Math.Exp(-q*T)*BS.NormCDF(g);

            double A = BSPut;
            double B = (v0*r0 + theta*r1)*P02;
            double C = (v0*p0 + theta*p1)*P11;
            double D = (v0*q0 + theta*q1)*P21 + 0.5*Math.Pow(v0*p0 + theta*p1,2)*P22;

            double[] output = new double[4];
            output[0] = A;
            output[1] = B;
            output[2] = C;
            output[3] = D;

            return output;
        }
    }
}

