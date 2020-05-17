using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Explicit_PDE_Method_American
{
    class ExplicitPDE
    {
        public double[,] HestonExplicitPDENonUniformGrid(HParam param,double K,double r,double q,double[] S,double[] V,double[] T,string PutCall,string EuroAmer)
        {
            // Finite differences for the Heston PDE for a European Call
            // Uses uneven grid sizes
            // In 'T Hout and Foulon "ADI Finite Difference Schemes for Option Pricing
            // in the Heston Modelo with Correlation" Int J of Num Analysis and Modeling, 2010.
            // INPUTS
            //    params = 6x1 vector of Heston parameters
            //    K = Strike price
            //    r = risk free rate
            //    q = dividend yield
            //    S = vector for stock price grid
            //    V = vector for volatility grid
            //    T = vector for maturity grid
            // OUTPUTS
            //    U = U(S,v) 2-D array of size (nS+1)x(nV+1)x(nT+1) for the call price

            // Heston parameters
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            // Grid measurements
            int NS = S.Length;
            int NV = V.Length;
            int NT = T.Length;
            double Smin = S[0]; double Smax = S[NS-1];

            double Vmin = V[0]; double Vmax = V[NV-1];
            double Tmin = T[0]; double Tmax = T[NT-1];
            double dt = (Tmax-Tmin)/Convert.ToDouble(NT-1);

            // Initialize the 2-D grid with zeros
            double[,] U = new double[NS,NV];

            // Temporary grid for previous time steps
            double[,] u = new double[NS,NV];

            //// Solve the PDE
            // Round each value of U(S,v) at each step
            // Boundary condition for t=maturity
            for(int s=0;s<=NS-1;s++)
                for(int v=0;v<=NV-1;v++)
                {
                    if (PutCall == "C")
                        U[s,v] = Math.Max(S[s] - K,0.0);
                    else
                        U[s,v] = Math.Max(K - S[s],0.0);
                }

            double LHS,derV,derS,derSS,derVV,derSV,L;

            for(int t=0;t<=NT-2;t++)
            {
                // Boundary condition for Smin and Smax
                for(int v=0;v<=NV-2;v++)
                {
                    U[0,v] = 0.0;
                    if (PutCall == "C")
                        U[NS-1,v] = Math.Max(0.0,Smax - K);
                    else
                        U[NS-1,v] = Math.Max(0.0,K - Smax);
                }
                // Boundary condition for Vmax
                for(int s=0;s<=NS-1;s++)
                {
                    if (PutCall == "C")
                        U[s,NV-1] = Math.Max(0.0,S[s] - K);
                    else
                        U[s,NV-1] = Math.Max(0.0,K - S[s]);
                }

                // Update the temporary grid u(s,t) with the boundary conditions
                for(int s=0;s<=NS-1;s++)
                    for(int v=0;v<=NV-1;v++)
                        u[s,v] = U[s,v];

                // Boundary condition for Vmin.
                // Previous time step values are in the temporary grid u(s,t)
                for(int s=1;s<=NS-2;s++)
                {
                    derV = (u[s,1]   - u[s,0])   / (V[1]-V[0]);          // Forward difference
                    derS = (u[s+1,0] - u[s-1,0]) / (S[s+1]-S[s-1]);      // Central difference
                    LHS = -r*u[s,0] + (r-q)*S[s]*derS + kappa*theta*derV;
                    U[s,0] = LHS*dt + u[s,0];
                }

                // Update the temporary grid u(s,t) with the boundary conditions
                for(int s=0;s<=NS-1;s++)
                    for(int v=0;v<=NV-1;v++)
                        u[s,v] = U[s,v];

                // Interior points of the grid (non boundary).
                // Previous time step values are in the temporary grid u(s,t)
                for(int s=1;s<=NS-2;s++)
                {
                    for(int v=1;v<=NV-2;v++)
                    {
                        derS  =  (u[s+1,v] - u[s-1,v]) / (S[s+1]-S[s-1]);  // Central difference for dU/dS
                        derV  =  (u[s,v+1] - u[s,v-1]) / (V[v+1]-V[v-1]);  // Central difference for dU/dV
                        derSS = ((u[s+1,v] - u[s,v])   / (S[s+1]-S[s]) - (u[s,v] - u[s-1,v])/(S[s]-S[s-1])) / (S[s+1]-S[s]);  // d2U/dS2
                        derVV = ((u[s,v+1] - u[s,v])   / (V[v+1]-V[v]) - (u[s,v] - u[s,v-1])/(V[v]-V[v-1])) / (V[v+1]-V[v]);  // d2U/dV2
                        derSV =  (u[s+1,v+1] - u[s-1,v+1] - U[s+1,v-1] + U[s-1,v-1]) / (S[s+1]-S[s-1]) / (V[v+1]-V[v-1]);     // d2U/dSdV
                        L = 0.5*V[v]*S[s]*S[s]*derSS + rho*sigma*V[v]*S[s]*derSV 
            			  + 0.5*sigma*sigma*V[v]*derVV - r*u[s,v]
			              + (r-q)*S[s]*derS + kappa*(theta-V[v])*derV;
                        U[s,v] = L*dt + u[s,v];
                    }
                }
                if(EuroAmer == "A")
                {
                    for(int s=0;s<=NS-1;s++)
                    {
                        for(int v=0;v<=NV-1;v++)
                        {
                            if(PutCall == "C")
                                U[s,v] = Math.Max(U[s,v],S[s] - K);
                            else
                                U[s,v] = Math.Max(U[s,v],K - S[s]);
                        }
                    }
                }
            }
            return U;
        }
    }
}

