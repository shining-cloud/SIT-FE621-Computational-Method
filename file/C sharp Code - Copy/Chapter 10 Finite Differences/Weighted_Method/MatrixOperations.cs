using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Weighted_Method
{
    class MatrixOps
    {
        // Approximate inverse of a matrix by Neumann series
        public double[,] MInvApprox(double[,] A,int M)
        {
            int n = A.GetLength(0);
            double[,] I = CreateI(n);
            double[,] Ainv = I;
            for(int k=1;k<=M;k++)
                Ainv = MMAdd(Ainv,MPower(MMSubtract(I,A),k));
            return Ainv;
        }
        // Inverse of a matrix through LU decomposition
        public double[,] MInvLU(double[,] A)
        {
            LUstruct Mats;
            Mats = LU(A);
            double[,] L = Mats.LM;
            double[,] U = Mats.UM;
            double[,] Uinv = MatUpTriangleInv(U);
            double[,] Linv = MatLowTriangleInv(L);
            double[,] Ainv = MMMult(Uinv,Linv);
            return Ainv;
        }

        // Multiply a matrix by a vector
        public double[] MVMult(double[,] A,double[] B)
        {
            int N = A.GetLength(0);
            int K = A.GetLength(1);
            double[] C = new double[N];
            for(int j=0;j<=N-1;j++)
            {
                C[j] = 0.0;
                for(int c=0;c<=K-1;c++)
                    C[j] += A[j,c] * B[c];
            }
            return C;
        }

        // Multiply two matrices together
        // First matrix is (n x k), Second is (k x m)
        // Resulting matrix is (n x m)
        public double[,] MMMult(double[,] A,double[,] B)
        {
            int N = A.GetLength(0);
            int K = A.GetLength(1);
            int M = B.GetLength(1);
            double[,] C = new double[N,M];

            for(int j=0;j<=M-1;j++)
            {
                for(int i=0;i<=N-1;i++)
                {
                    C[i,j] = 0.0;
                    for(int r=0;r<=K-1;r++)
                        C[i,j] += A[i,r] * B[r,j];
                }
            }
            return C;
        }

        // Transpose of a (N x M ) matrix
        public double[,] MTrans(double[,] A)
        {
            int N = A.GetLength(0);
            int M = A.GetLength(1);
            double[,] C = new double[M,N];
            for(int j=0;j<=M-1;j++)
            {
                for(int i=0;i<=N-1;i++)
                    C[j,i] = A[i,j];
            }
            return C;
        }

        // Trace of a matrix
        public double MTrace(double[,] A)
        {
            double sum = 0.0;
            int N = A.GetLength(0);
            for(int j=0;j<=N-1;j++)
                sum += A[j,j];
            return sum;
        }

        // Mean of a vector
        public double VMean(double[] X)
        {
            double XBar = 0.0;
            int N = X.Length;
            for(int i=0;i<=N-1;i++)
                XBar += X[i];
            return XBar / Convert.ToDouble(N);
        }

        // LU decomposition;        
        static LUstruct LU(double[,] A)
        {
            int N = A.GetLength(0);
            double[,] B = new double[N,N];
            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=N-1;j++)
                    B[i,j] = A[i,j];

            for(int k=0;k<=N-2;k++)
            {
                for(int i=k+1;i<=N-1;i++)
                    B[i,k] = B[i,k] / B[k,k];
                for(int j=k+1;j<=N-1;j++)
                {
                    for(int i=k+1;i<=N-1;i++)
                        B[i,j] = B[i,j] - B[i,k]*B[k,j];
                }
            }
            double[,] L = new double[N,N];
            double[,] U = new double[N,N];
            for(int i=0;i<=N-1;i++)
            {
                L[i,i] = 1.0;
                for(int j=0;j<=N-1;j++)
                {
                    if(i>j)
                        L[i,j] = B[i,j];
                    else
                        U[i,j] = B[i,j];
                }
            }
            LUstruct Mats;
            Mats.LM = L;
            Mats.UM = U;
            return Mats;
        }
        // Inverse of an upper triangular matrix
        public double[,] MatUpTriangleInv(double[,] U)
        {
            int N = U.GetLength(0);
            double[,] V = new double[N,N];
            for(int j=N-1;j>=0;j--)
            {
                V[j,j] = 1.0/U[j,j];
                for(int i=j-1;i>=0;i--)
                    for(int k=i+1;k<=j;k++)
                        V[i,j] -= 1.0 / U[i,i] * U[i,k] * V[k,j];

            }
            return V;
        }
        // Inverse of a lower triangular matrix
        public double[,] MatLowTriangleInv(double[,] L)
        {
            int N = L.GetLength(0);
            double[,] V = new double[N,N];
            for(int i=0;i<=N-1;i++)
            {
                V[i,i] = 1.0/L[i,i];
                for(int j=i-1;j>=0;j--)
                    for(int k=i-1;k>=j;k--)
                        V[i,j] -= 1.0 / L[i,i] * L[i,k] * V[k,j];
            }
            return V;
        }
        // Create the identity matrix
        public double[,] CreateI(int n)
        {
            double[,] I = new double[n,n];
            for(int i=0;i<=n-1;i++)
                for(int j=0;j<=n-1;j++)
                    if(i==j)
                        I[i,j] = 1.0;
                    else
                        I[i,j] = 0.0;
            return I;
        }
        // Power of a matrix A to the power k
        public double[,] MPower(double[,] A,int k)
        {
            int n = A.GetLength(0);
            double[,] Ak = new double[n,n];
            if(k==0)
                // Identity matrix
                Ak = CreateI(n);
            else
                // Recursive relation.  Faster to hold this off to higher powers
                Ak = MMMult(A,MPower(A,k-1));
            return Ak;
        }
        // Add two matrices each of size (nxk) together
        public double[,] MMAdd(double[,] A,double[,] B)
        {
            int N = A.GetLength(0);
            int K = A.GetLength(1);
            double[,] C = new double[N,K];
            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=K-1;j++)
                    C[i,j] = A[i,j] + B[i,j];
            return C;
        }
        // Substract two matrices each of size (nxk) together
        public double[,] MMSubtract(double[,] A,double[,] B)
        {
            int N = A.GetLength(0);
            int K = A.GetLength(1);
            double[,] C = new double[N,K];
            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=K-1;j++)
                    C[i,j] = A[i,j] - B[i,j];
            return C;
        }
        // Infinity Matrix Norm
        public double MNorm(double[,] A)
        {
            int N = A.GetLength(0);
            double[] rows = new double[N];
            for(int i=0;i<=N-1;i++)
            {
                rows[i] = 0.0;
                for(int j=0;j<=N-1;j++)
                    rows[i] += Math.Abs(A[i,j]);
            }
            return rows.Max();
        }
    }
}



