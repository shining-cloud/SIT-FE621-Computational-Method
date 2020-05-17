using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Moment_Matching
{
    class RandomNumbers
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
        // Box Muller transformation for simulating standard normal random variables
        public double RandomNorm()
        {
            double U1 = RandomNum(0.0,1.0);
            double U2 = RandomNum(0.0,1.0);
            return Math.Sqrt(-2.0*Math.Log(U1)) * Math.Sin(2.0*Math.PI*U2);
        }
        // Mean
        public double Mean(double[] x)
        {
            int M = x.Length;
            double N = Convert.ToDouble(M);
            return x.Sum()/N;
        }
    }
}

