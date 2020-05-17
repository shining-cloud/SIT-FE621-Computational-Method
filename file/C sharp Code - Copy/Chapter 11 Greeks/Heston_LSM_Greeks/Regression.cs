using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_LSM_Greeks
{
    class Regression
    {
        // Regression parameters
        public double[] Beta(double[,] X,double[] y)
        {
            double[,] Xt = MTrans(X);
            double[,] XtX1 = MInv(MMMult(Xt,X));
            double[] Xty   = MVMult(Xt,y);
            return MVMult(XtX1,Xty);
        }

        // Multiply matrix by a vector
        // Matrix is (n x k), vector is (k x 1)
        // Resulting matrix is (n x m)
        public double[] MVMult(double[,] A,double[] B)
        {
            int N = A.GetLength(0);
            int K = A.GetLength(1);
            double[] C = new double[N];
            for(int j=0;j<=N-1;j++)
            {
                C[j] = 0;
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
                    C[i,j] = 0;
                    for(int r=0;r<=K-1;r++)
                        C[i,j] += A[i,r] * B[r,j];
                }
            }
            return C;
        }

        // Minor of matrix A, eliminating row RowExclude and column ColExclude
        // RowExclude and ColExclude run from 0 to (n-1)
        public double[,] MMinor(double[,] A,int RowExclude,int ColExclude)
        {
            int N = A.GetLength(0);
            double[,] Minor = new double[N-1,N-1];
            int ColCount=0;
            for(int j=0;j<=N-1;j++)
            {
                if(j != ColExclude)
                {
                    int RowCount = 0;
                    for(int i=0;i<=N-1;i++)
                    {
                        if(i != RowExclude)
                        {
                            Minor[RowCount,ColCount] = A[i,j];
                            RowCount++;
                        }
                    }
                    ColCount++;
                }
            }
            return Minor;
        }

        // Determinant of a square matrix
        // Expand along the first row
        public double MDet(double[,] A)
        {
            int N = A.GetLength(0);
            double det = 0;
            if(N == 1) det = A[0,0];
            else
            {
                if(N == 2) det = A[0,0]*A[1,1] - A[1,0]*A[0,1];
                else
                {
                    for(int j=0;j<=N-1;j++)
                        det += A[0,j]*Math.Pow(-1,j)*MDet(MMinor(A,0,j));
                }
            }
            return det;
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

        // Inverse of a Matrix
        public double[,] MInv(double[,] A)
        {
            int N = A.GetLength(0);
            double[,] C = new double[N,N];
            double det = MDet(A);
            for(int j=0;j<=N-1;j++)
            {
                for(int i=0;i<=N-1;i++)
                    C[i,j] = Math.Pow(-1,(i+j))*MDet(MMinor(A,i,j)) / det;
            }
            return C;
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
    }
}

