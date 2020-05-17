using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Bounds_on_Alpha
{
    class AlphaBounds
    {
        static void Main(string[] args)
        {
            double S = 1.0;	    			    // Spot Price
            double K = 1.2; 				    // Strike Price
            double T = 1.0; 			        // Maturity in Years
            double rf = 0.0; 					// Interest Rate
            double q = 0.0;                     // Dividend yield
            double rho = -0.7;					// Heston Parameter: Correlation
            double kappa = 1.0;					// Heston Parameter 
            double theta = 0.1; 				// Heston Parameter 
            double lambda = 0.0;				// Heston Parameter 
            double sigma = 1.0;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.1;    				// Heston Parameter: Current Variance
            int trap = 1;                       // 1="Little Trap" characteristic function

            // Verify : kappa - rho*sigma >0
            double check = kappa - rho*sigma;
            if(check < 0)
            {
                Console.WriteLine("Warning: kappa - rho*sigma = {0:F3} is < 0",check);
                Console.WriteLine("------------------------------------------");
            }

            // Bounds on alpha and optimal alpha 
            // Shows multiple solutions to Rogers Lee's formula 
            // in which we solve for "a" in g(-ia)*exp(d(-ia)*T) = 1
            LordKahl LK = new LordKahl();
            double start = -15.0;
            double increment = 0.01;
            Complex g,d,E;
            double E1;
            Complex i = new Complex(0.0,1.0);
            double[] A = new Double[3001];
            double[] y = new Double[3001];
            double[] z = new Double[3001];
            for(int j=0;j<=3000;j++)
            {
                A[j] = start + j*increment;
                g = LK.RogerLeeG(-i*A[j],kappa,rho,sigma);
                d = LK.RogerLeeD(-i*A[j],kappa,rho,sigma);
                E = g*Complex.Exp(d*T);
                E1 = E.Real;
                y[j] = (E1 - 1.0);
                z[j] = 0.0;
            }
            // Find yMax and yMin from Roger Lee's closed form
            double ymax = (sigma - 2.0*kappa*rho + Math.Sqrt(sigma*sigma - 4.0*kappa*rho*sigma + 4.0*kappa*kappa))
	                    / (2.0*sigma*(1.0-rho*rho));
            double ymin = (sigma - 2.0*kappa*rho - Math.Sqrt(sigma*sigma - 4.0*kappa*rho*sigma + 4.0*kappa*kappa))
                        / (2.0*sigma*(1-rho*rho));

            // Settings for the Nelder Mead algorithm
            int N = 1;						// Number of parameters in f(x) to find
            int NumIters = 1;			    // First Iteration
            double MaxIters = 5000;			// Maximum number of iterations
            double Tolerance = 1e-20;		// Tolerance on best and worst function values

            // Starting values (vertices) in vector form.   Add more as needed
            double [,] s1 = new double [N,N+1];
            // Vertice 0	          Vertice 1
            s1[0,0] = ymin - 0.09;    s1[0,1] = ymin - 1.01;

            // Select the Roger Lee function as the objective function
            GlobalVars.ObjFunChoice = "RogerLee";

            // Arrange parameters in a vector to pass to Nelder Mead function
            double[] paramRL = new double[4];
            paramRL[0] = kappa;
            paramRL[1] = rho;
            paramRL[2] = sigma;
            paramRL[3] = T;

            // Calculate lower limit of the the range Ax
            NelderMeadAlgo NM = new NelderMeadAlgo();
            double[] AxLo = NM.NelderMead(NM.f,N,NumIters,MaxIters,Tolerance,s1,paramRL);
            double yneg = AxLo[0];

            // Calculate upper limit of the the range Ax
            double[,] s2 = new double[N,N+1];
            s2[0,0] = ymax + 4.99; s2[0,1] = ymax + 5.01;
            double[] AxHi = NM.NelderMead(NM.f,N,NumIters,MaxIters,Tolerance,s2,paramRL);
            double ypos = AxHi[0];

            // Bounds on alpha
            double AlphaMax = ypos - 1.0;
            double AlphaMin = yneg - 1.0;

            // Select the KahlLord function as the objective function
            GlobalVars.ObjFunChoice = "LordKahl";

            // Starting values for Lord and Kahl
            double StartVal = (AlphaMax + AlphaMin) / 2.0;
            double[,] s3 = new double[N,N+1];
            s3[0,0] = StartVal + 0.01;     s3[0,1] = StartVal - 0.01;

            // Parameters for Lord and Kahl
            double[] paramLK = new double[10];
            paramLK[0] = kappa;
            paramLK[1] = theta;
            paramLK[2] = lambda;
            paramLK[3] = rho;
            paramLK[4] = sigma;
            paramLK[5] = T;
            paramLK[6] = K;
            paramLK[7] = S;
            paramLK[8] = rf;
            paramLK[9] = v0;

            // Lord and Kahl optimal alpha
            double[] AlphaOptimal = NM.NelderMead(NM.f,N,NumIters,MaxIters,Tolerance,s3,paramLK);

            // Write the results
            Console.WriteLine("Roger Lee bounds on alpha");
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Ymin and Ymax are          ({0,7:F4}, {1,7:F4})",ymin,ymax);
            Console.WriteLine("The range Ax is            ({0,7:F4}, {1,7:F4})",yneg,ypos);
            Console.WriteLine("(AlphaMin,AlphMax) is      ({0,7:F4}, {1,7:F4}) = Ax - 1",AlphaMin,AlphaMax);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Lord and Kahl optimal Alpha {0,7:F4}", AlphaOptimal[0]);
            Console.WriteLine("--------------------------------------------------------");
            Console.WriteLine("Note that optimal alpha in ({0,7:F4}, {1,7:F4})",AlphaMin,AlphaMax);
            Console.WriteLine(" ");
        }
    }
}


