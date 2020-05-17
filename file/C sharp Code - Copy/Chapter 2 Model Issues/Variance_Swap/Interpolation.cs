using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Variance_Swap
{
    class Interpolate
    {
        // One dimensional linear interpolation
        public double interp1(double[] X,double[] Y,double xi)
        {
            int N = Y.Length;
            int x1 = 0;
            int x2 = 0;
            double yi = 0.0;

            // Look for xi on the end points
            if(xi == X[0])
                yi = Y[0];
            else if(xi == X[N-1])
                yi = Y[N-1];
            else
                for(int i=1;i<=N-1;i++)
                    if((X[i-1] <= xi) & (xi < X[i]))
                    {
                        x1 = i-1;
                        x2 = i;
                        double p = (xi - Convert.ToDouble(X[x1])) / (Convert.ToDouble(X[x2]) - Convert.ToDouble(X[x1]));
                        yi = (1-p)*Y[x1] + p*Y[x2];
                    }
            return yi;
        }
    }
}


