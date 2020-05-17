using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Differential_Evolution
{
    class NelderMeadAlgo
    {
        // Declare the form of the objective function ==========================================================================
        public delegate double ObjFun(double[] param,OFSet ofset);

        // Nelder Mead Algorithm ===============================================================================================
        public double[] NelderMead(ObjFun f,NMSet nmset,double[,] x)

        {
            int NumIters = 0;
            int N = nmset.N;
            int i,j;
            OFSet ofset = nmset.ofset;
            int MaxIters = nmset.MaxIters;
            double Tolerance = nmset.Tolerance;

            // Value of the function at the vertices
            double[][] F = new Double[N+1][];
            for(i=0;i<=N;i++)
                F[i] = new double[2] { 0.0,0.0 };

            // Step 0.  Ordering and Best and Worst points
        // Order according to the functional values, compute the best and worst points
        step0:
            NumIters = NumIters + 1;
            Console.Write("Nelder Mead iteration ");
            Console.WriteLine(NumIters);
            for(j=0;j<=N;j++)
            {
                double[] z = new double[N];
                for(i=0;i<=N-1;i++)
                {
                    z[i] = x[i,j];
                    F[j][0] = f(z,ofset);             // Function values
                    F[j][1] = j;											// Original index positions
                }
            }
            // Sort the F array w.r.t column 0
            int column = 0;
            Array.Sort(F,delegate(double[] w1,double[] w2)
            {
                return (w1[column] as IComparable).CompareTo(w2[column]);
            });

            // New vertices order first N best initial vectors and
            // last (N+1)st vertice is the worst vector
            // y is the matrix of vertices, ordered so that the worst vertice is last
            double[,] y = new double[N,N+1];
            for(j=0;j<=N;j++)
                for(i=0;i<=N-1;i++)
                    y[i,j] = x[i,Convert.ToInt32(F[j][1])];

            //  First best vector y(1) and function value f1
            double[] x1 = new double[N]; for(i=0;i<=N-1;i++) x1[i] = y[i,0];
            double f1 = f(x1,ofset);

            // Last best vector y(N) and function value fn
            double[] xn = new Double[N]; for(i=0;i<=N-1;i++) xn[i] = y[i,N-1];
            double fn = f(xn,ofset);

            // Worst vector y(N+1) and function value fn1
            double[] xn1 = new Double[N]; for(i=0;i<=N-1;i++) xn1[i] = y[i,N];
            double fn1 = f(xn1,ofset);

            // z is the first N vectors from y, excludes the worst y(N+1)
            double[,] zz = new Double[N,N];
            for(j=0;j<=N-1;j++)
                for(i=0;i<=N-1;i++) zz[i,j] = y[i,j];

            // Mean of best N values and function value fm
            double[] xm = new Double[N]; xm = VMean(zz,N);
            double fm = f(xm,ofset);

            // Reflection point xr and function fr
            double[] xr = new Double[N]; xr = VSub(VAdd(xm,xm),xn1);
            double fr = f(xr,ofset);

            // Expansion point xe and function fe
            double[] xe = new Double[N]; xe = VSub(VAdd(xr,xr),xm);
            double fe = f(xe,ofset);

            // Outside contraction point and function foc
            double[] xoc = new Double[N]; xoc = VAdd(VMult(xr,0.5),VMult(xm,0.5));
            double foc = f(xoc,ofset);

            // Inside contraction point and function foc
            double[] xic = new Double[N]; xic = VAdd(VMult(xm,0.5),VMult(xn1,0.5));
            double fic = f(xic,ofset);

            while((NumIters <= MaxIters) && (Math.Abs(f1-fn1) >= Tolerance))
            {
                // Step 1. Reflection Rule
                if((f1<=fr) && (fr<fn))
                {
                    for(j=0;j<=N-1;j++)
                        for(i=0;i<=N-1;i++) x[i,j] = y[i,j];
                    for(i=0;i<=N-1;i++) x[i,N] = xr[i];
                    goto step0;
                }
                // Step 2.  Expansion Rule
                if(fr<f1)
                {
                    for(j=0;j<=N-1;j++)
                        for(i=0;i<=N-1;i++) x[i,j] = y[i,j];
                    if(fe<fr)
                        for(i=0;i<=N-1;i++) x[i,N] = xe[i];
                    else
                        for(i=0;i<=N-1;i++) x[i,N] = xr[i];
                    goto step0;
                }
                // Step 3.  Outside contraction Rule
                if((fn<=fr) && (fr<fn1) && (foc<=fr))
                {
                    for(j=0;j<=N-1;j++)
                        for(i=0;i<=N-1;i++) x[i,j] = y[i,j];
                    for(i=0;i<=N-1;i++) x[i,N] = xoc[i];
                    goto step0;
                }
                // Step 4.  Inside contraction Rule
                if((fr>=fn1) && (fic<fn1))
                {
                    for(j=0;j<=N-1;j++)
                        for(i=0;i<=N-1;i++)
                            x[i,j] = y[i,j];
                    for(i=0;i<=N-1;i++) x[i,N] = xic[i];
                    goto step0;
                }
                // Step 5. Shrink Step
                for(i=0;i<=N-1;i++)
                    x[i,0] = y[i,0];

                for(i=0;i<=N-1;i++)
                    for(j=1;j<=N;j++)
                        x[i,j] = 0.5*(y[i,j] + x[i,0]);
                goto step0;
            }

            // Output component
            double[] outvec = new Double[N+2];
            for(i=0;i<=N-1;i++)
                outvec[i] = x1[i];
            outvec[N] = f1;
            outvec[N+1] = NumIters;
            return outvec;
        }

        // Vector functions ===================================================================================
        // Function to calculate the mean value of a set of N vectors each of dimension N
        // namely a (N x N) matrix
        static double[] VMean(double[,] X,int N)
        {
            double[] meanX = new Double[N];
            for(int i=0;i<=N-1;i++)
            {
                meanX[i]=0.0;
                for(int j=0;j<=N-1;j++)
                    meanX[i] += X[i,j];
                meanX[i] = meanX[i] / Convert.ToDouble(N);
            }
            return meanX;
        }
        // Function to add two vectors
        public double[] VAdd(double[] X,double[] Y)
        {
            int N = X.Length;
            double[] Z = new Double[N];
            for(int j=0;j<=N-1;j++)
                Z[j] = X[j] + Y[j];
            return Z;
        }
        // Function to subtract two vectors
        public double[] VSub(double[] X,double[] Y)
        {
            int N = X.Length;
            double[] Z = new Double[N];
            for(int j=0;j<=N-1;j++)
                Z[j] = X[j] - Y[j];
            return Z;
        }
        // Function to multiply a vector by a constant
        public double[] VMult(double[] X,double a)
        {
            int N = X.Length;
            double[] Z = new Double[N];
            for(int j=0;j<=N-1;j++)
                Z[j] = a*X[j];
            return Z;
        }
    }
}
