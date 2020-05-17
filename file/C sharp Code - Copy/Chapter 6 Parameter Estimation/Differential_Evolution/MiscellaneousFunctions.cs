using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Differential_Evolution
{
    class MiscFunctions
    {
        // Random number in (a,b) ==========================================================
        public readonly Random U = new Random();
        public readonly object sync = new object();
        public double RandomNum(double a,double b)
        {
            int divisor = 1000000000;
            lock(sync) { return a + (b-a)*U.Next(0,divisor)/divisor; }
        }

        // Random integer in (a,b) ==========================================================
        public readonly Random U1 = new Random();
        public readonly object sync1 = new object();
        public int RandomInt(int a,int b)
        {
            lock(sync1) { return U1.Next(a,b); }
        }

        // Random permutation of a vector of integers  =======================================
        public int[] RandomPerm(int[] a)
        {
            int N = a.Length;
            int[][] F = new int[N][];
            for(int i=0;i<=N-1;i++)
                F[i] = new int[2] { 0,0 };
            for(int j=0;j<=N-1;j++)
            {
                for(int i=0;i<=N-1;i++)
                {
                    F[j][0] = RandomInt(0,100);
                    F[j][1] = j;
                }
            }
            // Sort the F array w.r.t column 0
            int column = 0;
            Array.Sort(F,delegate(int[] w1,int[] w2)
            {
                return (w1[column] as IComparable).CompareTo(w2[column]);
            });
            int[] b = new int[N];
            for(int j=0;j<=N-1;j++)
                b[j] = F[j][1];
            return b;
        }

        // Vector of integers with one of the indices removed
        public int[] RemoveIndex(int[] a,int position)
        {
            int N = a.Length;
            int[] b = new int[N-1];
            if(position == 0)
            {
                for(int i=1;i<=N-1;i++)
                    b[i-1] = a[i];
            }
            else if(position == N-1)
            {
                for(int i=0;i<=N-2;i++)
                    b[i] = a[i];
            }
            else
            {
                for(int i=0;i<=position-1;i++)
                    b[i] = a[i];
                for(int i=position;i<=N-2;i++)
                    b[i] = a[i+1];
            }
            return b;
        }
    }
}

