using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GenerateGaussLobatto
{
    class GenerateGaussLobatto
    {
        static void Main(string[] args)
        {
            // Choose the number of points and generate the abscissas and weights
            int n = 32;

            // Tolerance and max number of iterations for the bisection algorithm
            double Tol = 1e-10;
            int MaxIter = 5000;

            // Starting and ending values for intervals containing roots, and # points
            double a = -1.0;
            double b = 1.0;
            int nI = 1500;

            // Gauss Laguerre abscissas and weights
            XW xw = GaussLobatto(n,a,b,nI,Tol,MaxIter);
            double[] x = xw.abscissas;
            double[] w = xw.weights;
            int nX = x.Length;

            // Output the results
            Console.WriteLine("{0:F0}-point Gauss Lobatto",x.Length);
            Console.WriteLine();
            Console.WriteLine(" Number  Abscissa   Weight");
            Console.WriteLine("------------------------------");
            for(int i=0;i<=nX-1;i++)
                Console.WriteLine("{0,5:F0} {1,10:F4} {2,10:F4}",i+1,x[i],w[i]);
            Console.WriteLine("------------------------------");
        }

        // Gauss Legendre abscissas and weights
        static XW GaussLobatto(int N,double a,double b,int nI,double Tol,int MaxIter)
        {
            // The Legendre polynomial
            int n = N-1;
            int m = Convert.ToInt16(Math.Floor(n/2.0));
            double[] L = new double[m+1];
            double[] dL = new double[m+1];
            for (int k=0; k<=m; k++)
            {
             	L[k]  = Math.Pow(0.5,n) * Math.Pow(-1.0,k) * factorial(2*n-2*k) / factorial(k)
 		              / factorial(n-k) / factorial(n-2*k);
	            dL[k] = Math.Pow(0.5,n) * Math.Pow(-1.0,k) * factorial(2*n-2*k) / factorial(k)
 		              / factorial(n-k) / factorial(n-2*k) * (n-2*k);
            }

            // Fill in the blank powers of L
            double[] P = new double[n+1];
            double[] dP = new double[n+1];
            for(int k=0;k<=n;k++)
            {
                if((k+1) % 2 == 0)
                {
                    P[n-k] = 0.0;
                    dP[n-k] = dL[k/2];
                }
                else
                {
                    P[n-k] = L[(k+1)/2];
                    dP[n-k] = 0.0;
                }
            }

            // Find the roots
            double[] x = findroot(dP,a,b,nI,Tol,MaxIter);

            // When n is odd, set the middle root to exactly zero
            if((n-1) % 2 == 1)
                x[(n-1)/2] = 0.0;

            // Create the abscissas by augmenting the endpoints with -1 and +1
            int R = x.Length;
            double[] X = new double[R+2];
            X[0] = -1.0;
            X[R+1] = 1.0;
            for(int k=1;k<=R;k++)
                X[k] = x[k-1];

            // Find the weights
            double[] Poly = new double[n+1];
            double[] w = new double[n+1];
            double sumPoly;
            for(int j=1;j<=n;j++) 
            {
                sumPoly = 0.0;
                for(int k=0;k<=n;k++)
                {
                    Poly[k] = P[n-k] * Math.Pow(X[j],n-k);
                    sumPoly += Poly[k];
                }
                w[j] = 2.0 / N / (N-1) / (sumPoly*sumPoly);
            }
            w[0] = 2.0 / N / (N-1);
            w[n] = w[0];

            // Output the results
            XW output;
            output.abscissas = X;
            output.weights = w;
            return output;
        }

        // Find roots of a polynomial
        static double[] findroot(double[] C,double a,double b,int nI,double Tol,int MaxIter)
        {
            // Find the Sturm sequence
            double[][] Sturm = sturm(C);
            int nS = C.Length;

            // Find the signs of the Sturm sequences over intervals
            SturmRoots sr = findintervals(a,b,nI,C,Sturm,nS);

            // Number of roots and start/end of the intervals
            int NRoots = sr.NRoots;
            List<double> StartInt = sr.StartInt;
            List<double> EndInt = sr.EndInt;

            // Apply the bisection algorithm to find the roots
            double[] root = new double[NRoots];
            for(int i=0;i<=NRoots-1;i++)
            {
                a = StartInt[i];
                b = EndInt[i];
                root[i] = Bisection(C,a,b,Tol,MaxIter);
            }
            return root;
        }
        // Find the intervals containing the roots of the polynomials using Sturm sequences
        static SturmRoots findintervals(double a,double b,int nI,double[] C,double[][] Sturm,int nS)
        {
            // Define the intervals along the x-axis
            double[] I = new double[nI];

            double increment = (b - a) / Convert.ToDouble(nI);
            for(int i=0;i<=nI-1;i++)
                I[i] = a + Convert.ToDouble(i)*increment;

            // Find the sign changes of the Sturm polynomials
            int M,ListSize = 0,Roots = 0;
            int[] Count = new int[nI];
            int[,] SignCoeff = new int[nI,nS];
            double coeff = 0.0;
            for(int i=0;i<=nI-1;i++)
            {
                for(int j=0;j<=nS-1;j++)
                {
                    coeff = 0.0;
                    M = Sturm[j].Length;
                    for(int k=0;k<=M-1;k++)
                    {
                        coeff += Sturm[j][k] * Math.Pow(I[i],Convert.ToDouble(k));
                    }
                    // Sign of the Sturm sequence (-1, 0, or +1)
                    SignCoeff[i,j] = Math.Sign(coeff);
                }
                // Remove zero signs from the Sturm sequence
                List<int> Signs = new List<int>();
                ListSize = 0;
                for(int j=0;j<=nS-1;j++)
                {
                    // Find the size of each sign list without zeros
                    if(SignCoeff[i,j] != 0)
                    {
                        Signs.Add(SignCoeff[i,j]);
                        ListSize += 1;
                    }
                }
                for(int k=1;k<=ListSize-1;k++)
                {
                    // Count the number of sign changes
                    if(Signs[k] != Signs[k-1])
                        Count[i] += 1;
                }
                Signs.Clear();
            }
            // Find the number of roots and intervals containing the roots
            List<double> StartInt = new List<double>();
            List<double> EndInt = new List<double>();
            for(int i=1;i<=nI-1;i++)
                if(Count[i-1] != Count[i])
                {
                    Roots += 1;
                    StartInt.Add(I[i-1]);
                    EndInt.Add(I[i]);
                }

            // Output the results into the structure SturmRoots
            SturmRoots output;
            output.NRoots = Roots;
            output.StartInt = StartInt;
            output.EndInt = EndInt;
            return output;
        }
        // Remainder of dividing polynomial P by polynomial Q
        static double[] polyrem(double[] P,double[] Q)
        {
            // rem = remainder, quo = quotient
            int nP = P.Length - 1;
            int nQ = Q.Length - 1;
            while(nQ >= 0 && Q[nQ] == 0)
                nQ--;
            double[] rem = new double[nP+1];
            Array.Copy(P,rem,nP+1);
            double[] quo = new double[P.Length];
            for(int k=nP-nQ;k>=0;k--)
            {
                quo[k] = rem[nQ+k]/Q[nQ];
                for(int j=nQ+k-1;j>=k;j--)
                    rem[j] -= quo[k]*Q[j-k];
            }
            for(int j=nQ;j<=nP;j++)
                rem[j] = 0.0;
            return rem;
        }
        // Differentiate a polynomial
        static double[] polydiff(double[] C)
        {
            int N = C.Length;
            double[] dC = new double[N-1];
            for(int k=1;k<=N-1;k++)
                dC[k-1] = C[k]*k;
            return dC;
        }
        // Sturm sequence of polynomials
        static double[][] sturm(double[] p)
        {
            int Order = p.Length;
            double[][] P = new double[Order][];
            P[0] = p;
            P[1] = polydiff(p);
            for(int j=2;j<=Order-1;j++)
            {
                P[j] = polyrem(P[j-2],P[j-1]);
                for(int k=0;k<P[j].Length;k++)
                    P[j][k] = -P[j][k];
            }
            return P;
        }
        // Evaluate a polynomial
        static double polyeval(double[] C,double x)
        {
            int N = C.Length;
            double P = C[0];
            for(int i=1;i<=N-1;i++)
                P += C[i] * Math.Pow(x,Convert.ToDouble(i));
            return P;
        }
        // Factorial
        static double factorial(int n)
        {
            if(n==0)
                return 1.0;
            else
                return n*factorial(n-1);
        }
        // N choose k
        static double nchoosek(int n,int k)
        {
            double nck = 1;
            for(int i=1;i<=k;i++)
                nck *= Convert.ToDouble(n-k+i) / Convert.ToDouble(i);
            return nck;
        }
        // Bisection Algorithm
        static double Bisection(double[] C,double a,double b,double Tol,int MaxIter)
        {
            double c = 0.0;
            double fc = 0.0;
            double fa = polyeval(C,a);
            double fb = polyeval(C,b);
            for(int x=0;x<=MaxIter;x++)
            {
                c = (a+b)/2.0;
                fc = polyeval(C,c);
                if(Math.Abs(fc)<Tol)
                    break;
                else
                {
                    if(fa*fc < 0)
                        b = c;
                    else
                        a = c;
                }
            }
            return c;
        }
    }
}

struct SturmRoots
{
    public int NRoots;
    public List<double> StartInt;
    public List<double> EndInt;
}
struct XW
{
    public double[] abscissas;
    public double[] weights;
}

