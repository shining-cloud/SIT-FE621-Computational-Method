using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bivariate_Trapezoidal_Rule
{
    class MainProgram
    {
        static void Main(string[] args)
        {
            // Points at which to evaluate the standard normal bivariate CDF
            double x = -0.515;
            double y =  0.243;
            // Exact value of from Matlab
            double TrueValue = 0.180751981501178;

            // Number of integration points for the (X,Y) grid 
            int Nx = 1000;
            int Ny = 1000;

            // Lower integration points and integration grid
            double XLow = -5.0;
            double YLow = -5.0;
            double hx = (x - XLow)/Nx;
            double hy = (y - YLow)/Ny;
            double[] X = new double[Nx+1];
            double[] Y = new double[Ny+1];
            for(int j=0;j<=Nx;j++)
                X[j] = XLow + j*hx;
            for(int j=0;j<=Ny;j++)
                Y[j] = YLow + j*hy;

            // The CDF using the double trapezoidal rule, and the error
            double TrapValue = DoubleTrapz(X,Y);
            double error = TrueValue - TrapValue;

            // Output the results
            Console.WriteLine("Double Trapezoidal rule using (Nx,Ny) = ({0:F0},{1:F0}) points",Nx,Ny);
            Console.WriteLine("------------------------------------------------------------");
            Console.WriteLine("Value with Matlab    {0,15:F12}",TrueValue);
            Console.WriteLine("Value with C#        {0,15:F12}",TrapValue);
            Console.WriteLine("Error between both   {0,15:F12}",error);
            Console.WriteLine("X-value              {0,6:F3}",X[Nx]);
            Console.WriteLine("Y-value              {0,6:F3}",Y[Ny]);
            Console.WriteLine("------------------------------------------------------------");
        }
        // The double trapezoidal rule
        static double DoubleTrapz(double[] X,double[] Y)
        {
            int nX = X.Length;
            int nY = Y.Length;
            double a,b,c,d;
            double sumInt = 0.0;
            for(int y=1;y<=nY-1;y++)
            {
                a = Y[y-1];
                b = Y[y];
                for(int x=1;x<=nX-1;x++)
                {
                    c = X[x-1];
                    d = X[x];
                    double term1 = f(a,c) + f(a,d) + f(b,c) + f(b,d);
                    double term2 = f((a+b)/2.0,c) + f((a+b)/2.0,d) + f(a,(c+d)/2.0) + f(b,(c+d)/2.0);
                    double term3 = f((a+b)/2.0,(c+d)/2.0);
                    sumInt  += (b-a)*(d-c)/16.0*(term1 + 2.0*term2 + 4.0*term3);
                }
            }
            return sumInt;
        }
        // The standard normal bivariate CDF
        static double f(double x,double y)
        {
            double pi = Math.PI;
            return Math.Exp(-0.5*(x*x + y*y)) / 2.0 / pi;
        }
    }
}
