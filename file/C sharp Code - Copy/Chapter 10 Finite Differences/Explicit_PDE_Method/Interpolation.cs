using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;

namespace Explicit_PDE_Method
{
    class Interpolation
    {
        // Two dimensional linear interpolation
        public double interp2(double[] X,double[] Y,double[,] Z,double xi,double yi)
        {
            int NX = X.Length;
            int NY = Y.Length;
            int xflag = 0;
            int yflag = 0;
            int x1 = 0;
            int x2 = 0;
            int y1 = 0;
            int y2 = 0;
            double IntValue = 999.0;

            // Find the index for X
            if(xi == X[0])
            {
                x1 = 0;
                xflag = 1;
            }
            else if(xi == X[NX-1])
            {
                x1 = NX-1;
                xflag = 1;
            }
            else
                for(int i=1;i<=NX-1;i++)
                    if((X[i-1] <= xi) & (xi < X[i]))
                    {
                        x1 = i-1;
                        x2 = i;
                    }

            // Find the index for Y
            if(yi == Y[0])
            {
                y1 = 0;
                yflag = 1;
            }
            else if(yi == Y[NY-1])
            {
                y1 = NY-1;
                yflag = 1;
            }
            else
                for(int i=1;i<=NY-1;i++)
                    if((Y[i-1] <= yi) & (yi < Y[i]))
                    {
                        y1 = i-1;
                        y2 = i;
                    }

            // Interpolation: both xi and yi lie off the grid points
            if((xflag==0) & (yflag==0))
            {
                double z11 = Z[y1,x1];
                double z12 = Z[y1,x2];
                double z21 = Z[y2,x1];
                double z22 = Z[y2,x2];
                double px = (xi - Convert.ToDouble(X[x1])) / (Convert.ToDouble(X[x2]) - Convert.ToDouble(X[x1]));
                double py = (yi - Convert.ToDouble(Y[y1])) / (Convert.ToDouble(Y[y2]) - Convert.ToDouble(Y[y1]));
                double Y1int = (1.0-py)*z11 + py*z21;
                double Y2int = (1.0-py)*z12 + py*z22;
                IntValue = (1.0-px)*Y1int + px*Y2int;
            }
            // Interpolation: xi lies on the grid point, yi lies off the grid point
            if((xflag==1) &(yflag==0))
            {
                double z11 = Z[y1,x1];
                double z21 = Z[y2,x1];
                double py = (yi - Convert.ToDouble(Y[y1])) / (Convert.ToDouble(Y[y2]) - Convert.ToDouble(Y[y1]));
                IntValue = (1.0-py)*z11 + py*z21;
            }

            // Interpolation: xi lies off the grid point, yi lies on the grid point
            if((xflag==0) & (yflag==1))
            {
                double z11 = Z[y1,x1];
                double z12 = Z[y1,x2];
                double px = (xi - Convert.ToDouble(X[x1])) / (Convert.ToDouble(X[x2]) - Convert.ToDouble(X[x1]));
                IntValue = (1.0-px)*z11 + px*z12;
            }

            // Interpolation: both xi and yi lie on the grid;
            if((xflag==1) & (yflag==1))
                IntValue = Z[y1,x1];

            // Return the result
            return IntValue;
        }

        // Inverse hyperbolic sine
        public double aSinh(double x)
        {
            return Math.Log(x + Math.Sqrt(x*x+1.0));
        }
    }
}


