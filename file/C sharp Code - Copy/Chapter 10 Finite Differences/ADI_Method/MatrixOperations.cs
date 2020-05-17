using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ADI_Method
{
    class MatrixOps
    {
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

        // Inverse of a Matrix using determinants
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
            double det = 0.0;
            if(N == 1)
                det = A[0,0];
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

        // Multiply a matrix by a scalar
        public double[,] MMultS(double[,] A,double S)
        {
            int N = A.GetLength(0);
            int M = A.GetLength(1);
            double[,] B = new double[N,M];
            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=M-1;j++)
                    B[i,j] = S * A[i,j];
            return B;
        }
        // Multiply a vector by a scalar
        public double[] VMultS(double[] A,double S)
        {
            int N = A.Length;
            double[] B = new double[N];
            for(int i=0;i<=N-1;i++)
                    B[i] = S * A[i];
            return B;
        }

        // Add two matrices
        public double[,] MAdd(double[,] A,double[,] B)
        {
            int N = A.GetLength(0);
            int M = A.GetLength(1);
            double[,] C = new double[N,M];
            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=M-1;j++)
                    C[i,j] = A[i,j] + B[i,j];
            return C;
        }

        // Substract two matrices
        public double[,] MSub(double[,] A,double[,] B)
        {
            int N = A.GetLength(0);
            int M = A.GetLength(1);
            double[,] C = new double[N,M];
            for(int i=0;i<=N-1;i++)
                for(int j=0;j<=M-1;j++)
                    C[i,j] = A[i,j] - B[i,j];
            return C;
        }

        // Add two vectors
        public double[] VAdd(double[] A,double[] B)
        {
            int N = A.Length;
            double[] C = new double[N];
            for(int i=0;i<=N-1;i++)
                C[i] = A[i] + B[i];
            return C;
        }


        // Substract two vectors
        public double[] VSub(double[] A,double[] B)
        {
            int N = A.Length;
            double[] C = new double[N];
            for(int i=0;i<=N-1;i++)
                    C[i] = A[i] - B[i];
            return C;
        }
    }
}



