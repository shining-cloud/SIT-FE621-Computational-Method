using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Explicit_PDE_Method
{
    class PDEAlgoU
    {
        public double[,] HestonExplicitPDE(HParam param,double K,double r,double q,double[] S,double[] V,double[] T)
        {
            // Finite differences for the Heston PDE for a European Call
            // Uses even grid sizes
            // In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
            // in the Heston Modelo with Correlation" Int J of Num Analysis and Modeling, 2010.
            // Thesis by Sensi Li and paper by Vassilis Galiotos
            // INPUTS
            //    params = 6x1 vector of Heston parameters
            //    K = Strike price
            //    r = risk free rate
            //    q = dividend yield
            //    S = vector for stock price grid
            //    V = vector for volatility grid
            //    T = vector for maturity grid
            // OUTPUT
            //    U = U(S,v) 2-D array of size (nS+1)x(nV+1) for the call price

            // Heston parameters
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            // Length of stock price, volatility, and maturity
            int NS = S.Length;
            int NV = V.Length;
            int NT = T.Length;
            double Smin = S[0]; double Smax = S[NS-1];
            double Vmin = V[0]; double Vmax = V[NV-1];
            double Tmin = T[0]; double Tmax = T[NT-1];

            // Increment for Stock Price, Volatility, and Maturity
            double ds = (Smax-Smin)/Convert.ToDouble(NS-1);
            double dv = (Vmax-Vmin)/Convert.ToDouble(NV-1);
            double dt = (Tmax-Tmin)/Convert.ToDouble(NT-1);

            // Initialize the 2-D grid with zeros
            double[,] U = new double[NS,NV];

            // Temporary grid for previous time steps
            double[,] us = new double[NS,NV];

            double DerV,DerS,A,B,C,D,E,F;

            //// Solve the PDE
            // Round each value of U(S,v,t) at each step
            // Boundary condition for t = Maturity
            for(int s=0;s<=NS-1;s++)
                for(int v=0;v<=NV-1;v++)
                    U[s,v] = Math.Max(S[s] - K,0.0);

            for(int t=0;t<=NT-2;t++)
            {
                // Boundary condition for Smin and Smax
                for(int v=0;v<=NV-2;v++)
                {
                    U[0,v] = 0.0;
                    U[NS-1,v] = Math.Max(0.0,Smax - K); // Galiotos uses U(NS-1,v) + ds;
                }

                // Boundary condition for Vmax
                for(int s=0;s<=NS-1;s++)
                    U[s,NV-1] = Math.Max(0.0,S[s] - K);     // Galiotos uses U(s,NV-1);

                // Update the temporary grid U(s,t) with the boundary conditions
                for (int s=0;s<=NS-1; s++)
                    for (int v=0; v<=NV-1; v++)
                        us[s,v] = U[s,v];

                // Boundary condition for Vmin.
                // Previous time step values are in the temporary grid u(s,t)
                for(int s=1;s<=NS-2;s++)
                {
                    DerV = (us[s,1] - us[s,0]) / dv;                	        // PDE Points on the middle of the grid (non boundary)
                    DerS = (us[s+1,0] - us[s-1,0])/2.0/ds;                    // Central difference for dU/dS
                    U[s,0] = us[s,0]*(1.0 - r*dt - kappa*theta*dt/dv)
                           + dt*(0.5*(r-q)*Convert.ToDouble(s+1)*us[s+1,0] - us[s-1,0])
                           + kappa*theta*dt/dv*us[s,1];
                }
                // Update the temporary grid us(s,t) with the boundary conditions
                for(int s=0;s<=NS-1;s++)
                    for(int v=0;v<=NV-1;v++)
                        us[s,v] = U[s,v];

                // Interior points of the grid (non boundary).
                // Previous time step values are in the temporary grid us(s,t)
                for(int s=1;s<=NS-2;s++)
                    for(int v=1;v<=NV-2;v++)
                    {
                        double s1 = Convert.ToDouble(s);
                        double v1 = Convert.ToDouble(v);
                        A = (1.0 - dt*s1*s1*v1*dv - sigma*sigma*v1*dt/dv - r*dt);
                        B = (0.5*dt*s1*s1*v1*dv - 0.5*dt*(r-q)*s1);
                        C = (0.5*dt*s1*s1*v1*dv + 0.5*dt*(r-q)*s1);
                        D = (0.5*dt*sigma*sigma*v1/dv - 0.5*dt*kappa*(theta-v1*dv)/dv);
                        E = (0.5*dt*sigma*sigma*v1/dv + 0.5*dt*kappa*(theta-v1*dv)/dv);
                        F = 0.25*dt*sigma*s1*v1;
                        U[s,v] = A*us[s,v] + B*us[s-1,v] + C*us[s+1,v]
				               + D*us[s,v-1] + E*us[s,v+1]
				               + F*(us[s+1,v+1]+us[s-1,v-1]-us[s-1,v+1]-us[s+1,v-1]);
                    }
            }
            return U;
        }
    }
}


