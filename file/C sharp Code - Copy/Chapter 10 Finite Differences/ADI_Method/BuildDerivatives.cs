using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ADI_Method
{
    class BuildDerivativesU
    {
        public LMatrices BuildDerivatives(double[] S,double[] V,double[] T)
        {
            // Build the first- and second-order derivatives for the "L" operator matrix
            // for the Weighted method
            // Requires a uniform grid for S and V
            // INPUTS
            //    S = vector for uniform stock price grid
            //    V = vector for uniform volatility grid
            //    T = vector for uniform maturity grid
            //    thet = parameter for the weighted scheme
            // OUTPUTS
            //    Matrices of dimnension N x N (N=NS+NV)
            //    derS  = Matrix for first-order derivative dU/dS 
            //    derSS = Matrix for second-order derivative dU2/dS2 
            //    derV1 = Matrix for first-order derivative dU/dV, kappa*theta portion
            //    derV2 = Matrix for first-order derivative dU/dV, -kappa*V portion
            //    derVV = Matrix for first-order derivative dU2/dV2
            //    derSV = Matrix for first-order derivative dU2/dSdV
            //    R     = Matrix for r*S(v,t) portion of the PDE

            // Length of stock price, volatility, and maturity
            int NS = S.Length;
            int NV = V.Length;
            int NT = T.Length;
            double Smin = S[0]; double Smax = S[NS-1];
            double Vmin = V[0]; double Vmax = V[NV-1];
            double Tmin = T[0]; double Tmax = T[NT-1];

            // Increment for Stock Price, Volatility, and Maturity
            double ds = (Smax-Smin) / Convert.ToDouble(NS-1);
            double dv = (Vmax-Vmin) / Convert.ToDouble(NV-1);
            double dt = (Tmax-Tmin) / Convert.ToDouble(NT-1);

            // Preliminary quantities
            // Size of the U(t) vector and L matrix
            int N = NS*NV;

            // The vectors for S and V, stacked
            double[] Si = new double[N];
            double[] Vi = new double[N];
            int k = 0;
            for(int v=0;v<=NV-1;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    Si[k] = S[s];
                    Vi[k] = V[v];
                    k += 1;
                }

            // Identification of the boundary points
            int[] VminB = new int[N];
            int[] VmaxB = new int[N];
            int[] SminB = new int[N];
            int[] SmaxB = new int[N];

            // Vmin and Vmax
            for(int v=0;v<=NS-2;v++)
                VminB[v] = 1;
            for(int v=N-NS+1;v<=N-2;v++)
                VmaxB[v] = 1;
            VmaxB[N-1] = 1;

            // Smin and Smax
            k = 0;
            for(int v=0;v<=NV-1;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    if(s == 0)
                        SminB[k] = 1;
                    else if(s == NS-1)
                        SmaxB[k] = 1;
                    k += 1;
                }
            SminB[0] = 0;
            SmaxB[NS-1] = 0;
            SmaxB[N-1] = 0;

            // Identification of the non-boundary points
            int[] NB = new int[N];
            for(k=0;k<=N-1;k++)
                if(SminB[k]==0 & SmaxB[k]==0 & VminB[k]==0 & VmaxB[k]==0)
                    NB[k] = 1;

            // Forward, backward and central differences for S
            int[] Cs = new int[N];
            int[] Fs = new int[N];
            int[] Bs = new int[N];
            k = 0;
            for(int v=0;v<=NV-2;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    if(s == 1)
                        Fs[k] = 1;
                    else if(s == NS-2)
                        Bs[k] = 1;
                    else if(s>=2 & s<=NS-3)
                        Cs[k] = 1;
                    k += 1;
                }
            Fs[1] = 0;
            Bs[NS-2] = 0;
            for(k=2;k<=NS-3;k++)
                Cs[k] = 0;

            // Forward, backward and central differences for V
            int[] Cv = new int[N];
            int[] Fv = new int[N];
            int[] Bv = new int[N];
            for(k=NS+1;k<=2*NS-2;k++)
                Fv[k] = 1;
            for(k=(NV-2)*NS+1;k<=(NV-1)*NS-2;k++)
                Bv[k] = 1;
            k = 0;
            for(int v=0;v<=NV-3;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    if(s>=1 & s<=NS-2)
                        Cv[k] = 1;
                    k += 1;
                }
            for(k=0;k<=2*NS-1;k++)
                Cv[k] = 0;

            // Central difference for SV-derivatives
            int[] Csv = new int[N];
            k = 0;
            for(int v=0;v<=NV-2;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    if(s>=1 & s<=NS-2)
                        Csv[k] = 1;
                    k += 1;
                }
            for(k=0;k<=NS-1;k++)
                Csv[k] = 0;

            // for(k=0;k<=N-1;k++)
            //      Console.WriteLine("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}",SminB[k],SmaxB[k],VminB[k],VmaxB[k],NB[k],Fs[k],Bs[k],Cs[k],Fv[k],Bv[k],Cv[k],Csv[k]);

            // Create the matrices for the derivatives
            double[,] derS  = new double[N,N];
            double[,] derSS = new double[N,N];
            double[,] derV1 = new double[N,N];
            double[,] derV2 = new double[N,N];
            double[,] derVV = new double[N,N];
            double[,] derSV = new double[N,N];
            double[,] R     = new double[N,N];

            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=N-1;j++)
                    if(i==j)
                        R[i,j] = 1.0;

            int M;
            int[] I;

            // Create the matrices for the derivatives
            // NON BOUNDARY POINTS ----------------------------------
            I = Find(Cs);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {                                                                     // Central differences  
                derS[I[k],I[k]-1]  =  -0.5/ds * Si[I[k]];                         // U(s-1,v)
                derS[I[k],I[k]]    =   0.0;                                       // U(s,v)
                derS[I[k],I[k]+1]  =   0.5/ds * Si[I[k]];                         // U(s+1,v)
                derSS[I[k],I[k]-1] =   1.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];    // U(s-1,v)
                derSS[I[k],I[k]]   =  -2.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];    // U(s,v)
                derSS[I[k],I[k]+1] =   1.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];    // U(s+1,v)
            }
            I = Find(Fs);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {                                                                   // Forward differences
                derS[I[k],I[k]]    = -1.5/ds  * Si[I[k]];                       // U(s,v)
                derS[I[k],I[k]+1]  =  2.0/ds  * Si[I[k]];                       // U(s+1,v)
                derS[I[k],I[k]+2]  = -1.0/ds  * Si[I[k]];                       // U(s+2,v)
                derSS[I[k],I[k]]   =  1.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];   // U(s,v)
                derSS[I[k],I[k]+1] = -2.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];   // U(s+1,v)
                derSS[I[k],I[k]+2] =  1.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];   // U(s+2,v)
            }
            I = Find(Bs);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {                                                                       // Backward differences
                derS[I[k],I[k]-2]  =  -0.5/ds * Si[I[k]];                           // U(s-2,v)
                derS[I[k],I[k]-1]  =  -2.0/ds * Si[I[k]];                           // U(s-1,v)
                derS[I[k],I[k]]    =   1.5/ds * Si[I[k]];                           // U(s,v)
                derSS[I[k],I[k]-2] =   1.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];      // U(s-2,v)
                derSS[I[k],I[k]-1] =  -2.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];      // U(s-1,v)
                derSS[I[k],I[k]]   =   1.0/ds/ds * Vi[I[k]]*Si[I[k]]*Si[I[k]];      // U(s,v)
            }

            // Create the matrix for V-derivatives
            I = Find(Cv);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {                                                       // Central differences
                derV1[I[k],I[k]-NS]  = -0.5/dv;                   // U(s,v-1)
                derV1[I[k],I[k]]     =  0.0;                      // U(s,v)
                derV1[I[k],I[k]+NS]  =  0.5/dv;                   // U(s,v+1)
                derV2[I[k],I[k]-NS]  = -0.5/dv * Vi[I[k]];        // U(s,v-1)
                derV2[I[k],I[k]]     =  0.0;                      // U(s,v)
                derV2[I[k],I[k]+NS]  =  0.5/dv * Vi[I[k]];        // U(s,v+1)
                derVV[I[k],I[k]-NS]  =  1.0/dv/dv * Vi[I[k]];        // U(s,v-1)
                derVV[I[k],I[k]]     = -2.0/dv/dv * Vi[I[k]];        // U(s,v)
                derVV[I[k],I[k]+NS]  =  1.0/dv/dv * Vi[I[k]];        // U(s,v+1)
            }
            I = Find(Fv);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {                                                       // Forward differences
                derV1[I[k],I[k]]       = -1.5/dv;                 // U(s,v)
                derV1[I[k],I[k]+NS]    =  2.0/dv;                 // U(s,v+1)
                derV1[I[k],I[k]+2*NS]  = -1.0/dv;                 // U(s,v+2)
                derV2[I[k],I[k]]       = -1.5/dv * Vi[I[k]];      // U(s,v)
                derV2[I[k],I[k]+NS]    =  2.0/dv * Vi[I[k]];      // U(s,v+1)
                derV2[I[k],I[k]+2*NS]  = -1.0/dv * Vi[I[k]];      // U(s,v+2)
                derVV[I[k],I[k]]       =  1.0/dv/dv * Vi[I[k]];    // U(s,v)
                derVV[I[k],I[k]+NS]    = -2.0/dv/dv * Vi[I[k]];    // U(s,v+1)
                derVV[I[k],I[k]+2*NS]  =  1.0/dv/dv * Vi[I[k]];    // U(s,v+2)
            }
            I = Find(Bv);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {                                                       // Backward differences
                derV1[I[k],I[k]-2*NS]  = -0.5/dv;                 // U(s,v-2)
                derV1[I[k],I[k]-NS]    = -2.0/dv;                 // U(s,v-1)
                derV1[I[k],I[k]]       =  1.5/dv;                 // U(s,v)
                derV2[I[k],I[k]-2*NS]  = -0.5/dv * Vi[I[k]];      // U(s,v-2)
                derV2[I[k],I[k]-NS]    = -2.0/dv * Vi[I[k]];      // U(s,v-1)
                derV2[I[k],I[k]]       =  1.5/dv * Vi[I[k]];      // U(s,v)
                derVV[I[k],I[k]-2*NS]  =  1.0/dv/dv * Vi[I[k]];      // U(s,v-2)
                derVV[I[k],I[k]-NS]    = -2.0/dv/dv * Vi[I[k]];      // U(s,v-1)
                derVV[I[k],I[k]]       =  1.0/dv/dv * Vi[I[k]];      // U(s,v)
            }
            // Create the matrix for SV-derivatives
            I = Find(Csv);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {
                derSV[I[k],I[k]+NS+1] =  1.0/(4.0*ds*dv) * Vi[I[k]]*Si[I[k]];  // U(s+1,v+1)
                derSV[I[k],I[k]+NS-1] = -1.0/(4.0*ds*dv) * Vi[I[k]]*Si[I[k]];  // U(s-1,v+1)
                derSV[I[k],I[k]-NS-1] =  1.0/(4.0*ds*dv) * Vi[I[k]]*Si[I[k]];  // U(s-1,v-1)
                derSV[I[k],I[k]-NS+1] = -1.0/(4.0*ds*dv) * Vi[I[k]]*Si[I[k]];  // U(s+1,v-1)
            }

            // BOUNDARY POINTS ----------------------------------
            // Boundary for Smin
            I = Find(SminB);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {
                derS[I[k],I[k]]  = 0.0;
                derSS[I[k],I[k]] = 0.0;
                derV1[I[k],I[k]] = 0.0;
                derV2[I[k],I[k]] = 0.0;
                derVV[I[k],I[k]] = 0.0;
                derSV[I[k],I[k]] = 0.0;
                R[I[k],I[k]] = 0.0;
            }
            // Boundary condition for Smax
            I = Find(SmaxB);
            M = I.Length;
            for(k=0;k<=M-1;k++)
            {
                derS[I[k],I[k]]     = Si[I[k]];
                derSS[I[k],I[k]]    = 0.0;
                derSV[I[k],I[k]]    = 0.0;                   // Central difference
                derV1[I[k],I[k]-NS] = -0.5/dv;               // U(s,v-1)
                derV1[I[k],I[k]]    =  0;                    // U(s,v)
                derV1[I[k],I[k]+NS] =  0.5/dv;               // U(s,v+1)
                derV2[I[k],I[k]-NS] = -0.5/dv * Vi[I[k]];    // U(s,v-1)
                derV2[I[k],I[k]]    =  0.0;                  // U(s,v)
                derV2[I[k],I[k]+NS] =  0.5/dv * Vi[I[k]];    // U(s,v+1)
                derVV[I[k],I[k]-NS] =  1.0/dv/dv * Vi[I[k]];    // U(s,v-1)
                derVV[I[k],I[k]]    = -2.0/dv/dv * Vi[I[k]];    // U(s,v)
                derVV[I[k],I[k]+NS] =  1.0/dv/dv * Vi[I[k]];    // U(s,v+1)
            }

            // Boundary condition for Vmax.  Only the LS submatrix is non-zero
            I = Find(VmaxB);
            M = I.Length;
            for(k=0;k<=M-1;k++)
                derS[I[k],I[k]] = Si[I[k]];

            // Boundary condition for Vmin
            // First row
            derS[0,0]     = -1.5/ds  * Si[0];     // U(s,v)
            derS[0,1]     =  2.0/ds  * Si[0];     // U(s+1,v)
            derS[0,2]     = -1.0/ds  * Si[0];     // U(s+2,v)
            derV1[0,0]    = -1.5/dv;              // U(s,v)
            derV1[0,NS]   =  2.0/dv;              // U(s,v+1)
            derV1[0,2*NS] = -1.0/dv;              // U(s,v+2)
            // Last row
            derS[N-2,N-4]       = -0.5/ds * Si[N-2];      // U(s-2,v)
            derS[N-2,N-3]       = -2.0/ds * Si[N-2];      // U(s-1,v)
            derS[N-2,N-2]       =  1.5/ds * Si[N-2];      // U(s,v)

            derV1[NS-2,NS-2+2*NS] = -0.5/dv * Vi[NS-2];      // U(s,v+2)
            derV1[NS-2,NS-2+NS]   =  2.0/dv * Vi[NS-2];      // U(s,v+1)
            derV1[NS-2,NS-2]      = -1.5/dv * Vi[NS-2];      // U(s,v)
            // Other rows
            for(int i=1;i<=NS-3;i++)
            {
                derS[i,i-1]     = -0.5/ds * Si[i];      // U(s-1,v)
                derS[i,i]       =  0.0;                 // U(s,v)
                derS[i,i+1]     =  0.5/ds * Si[i];      // U(s+1,v)
                derV1[i,i]      = -1.5/dv;              // U(s,v)
                derV1[i,i+NS]   =  2.0/dv;              // U(s,v+1)
                derV1[i,i+2*NS] = -1.0/dv;              // U(s,v+2)
            }

            // Output the sub-matrices
            LMatrices LMat;
            LMat.derS = derS;
            LMat.derSS = derSS;
            LMat.derV1 = derV1;
            LMat.derV2 = derV2;
            LMat.derVV = derVV;
            LMat.derSV = derSV;
            LMat.R = R;

            return LMat;
        }
        // Function that returns positions of vector entries with "1"
        public int[] Find(int[] X)
        {
            int[] G;
            int N = X.Length;
            List<int> F = new List<int>();
            for(int k=0;k<=N-1;k++)
                if(X[k] == 1)
                    F.Add(k);

            G = F.ToArray();
            return G;
        }
    }
}

