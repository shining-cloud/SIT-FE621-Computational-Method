using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace Differential_Evolution
{
    class DiffEvoAlgo
    {
        // Differential evolution algorithm
        public HParam HestonDE(DEParam DEsettings,OPSet settings,MktData data,int LossFunction,double a,double b,double Tol,int MaxIter,double[] X,double[] W,string CF)
        {
            // NG = Number of generations (iterations)
            // NP = Number of population members
            // CR = Crossover ratio (=0.5)
            // F  = Threshold (=0.8)
            // Hi = Vector of upper bounds for the parameters
            // Lo = Vector of lower bounds for the parameters
            // S  = Spot Price
            // K1 = First strike point for the FRFT
            // rf = Risk free rate
            // q  = Dividend Yield
            // MktPrice = Market quotes for prices
            // K = Vector of Strikes
            // T = Vector of Maturities
            // PutCall = Matrix of 'P'ut or 'C'all
            // MktIV = Market quotes for implied volatilities
            // LossFunction = Type of Objective Function 1 = MSE; 2 = RMSE; 3 = IVMSE; 4 = Christoffersen, Jacobs, Heston (2009)
            // Bisection method settings; a = Lower limit; b = upper limit; Tol = Tolerance; MaxIter = Max number of iterations
            // trap = 1 is "Little Trap" c.f., 0 is Heston c.f.
            // X, W = Gauss Laguerre abscissas and weights
            // CF = "Heston" or "Attari" characteristic function

            ObjectiveFunction OF = new ObjectiveFunction();
            MiscFunctions MF = new MiscFunctions();
            double[] Hi = DEsettings.ub;
            double[] Lo = DEsettings.lb;
            double CR   = DEsettings.CR;
            double kappaU = Hi[0]; double kappaL = Lo[0];
            double thetaU = Hi[1]; double thetaL = Lo[1];
            double sigmaU = Hi[2]; double sigmaL = Lo[2];
            double v0U    = Hi[3]; double v0L    = Lo[3];
            double rhoU   = Hi[4]; double rhoL   = Lo[4];

            // Create the structure for the objective function;
            OFSet ofset;
            ofset.opset = settings;
            ofset.data = data;
            ofset.X = X;
            ofset.W = W;
            ofset.LossFunction = LossFunction;
            ofset.lb = DEsettings.lb;
            ofset.ub = DEsettings.ub;
            ofset.CF = CF;

            // Step1.  Generate the population matrix of random parameters
            int NP = DEsettings.NP;
            int NG = DEsettings.NG;
            double F = DEsettings.F;
            double[,] P = new double[5,NP];
            for(int j=0;j<=NP-1;j++)
            {
                P[0,j] = kappaL + (kappaU-kappaL)*MF.RandomNum(0.0,1.0);
                P[1,j] = thetaL + (thetaU-thetaL)*MF.RandomNum(0.0,1.0);
                P[2,j] = sigmaL + (sigmaU-sigmaL)*MF.RandomNum(0.0,1.0);
                P[3,j] = v0L    + (v0U   -   v0L)*MF.RandomNum(0.0,1.0);
                P[4,j] = rhoL   + (rhoU  -  rhoL)*MF.RandomNum(0.0,1.0);
            }

            // Generate the random numbers outside the loop
            double[, ,] U = new double[5,NP,NG];
            for(int k=0;k<=NG-1;k++)
                for(int j=0;j<=NP-1;j++)
                    for(int i=0;i<=4;i++)
                        U[i,j,k] = MF.RandomNum(0.0,1.0);

            // Initialize the variables
            double[] Pr1 = new double[5];
            double[] Pr2 = new double[5];
            double[] Pr3 = new double[5];
            double[] P0  = new double[5];

            // Loop through the generations
            for(int k=0;k<=NG-1;k++)
            {
                Console.Write("Differential Evolution iteration "); Console.WriteLine(k);

                // Loop through the population
                for(int i=0;i<=NP-1;i++)
                {
                    // Select the i-th member of the population
                    for(int s=0;s<=4;s++)
                        P0[s] = P[s,i];

                    Step0:
                        // Select random indices for three other distinct members
                        int[] Integers0toNP1 = new int[NP];
                        for(int s=0;s<=NP-1;s++)
                            Integers0toNP1[s] = s;


                        // Random perumation of the indices (0,1,...,NP-1)
                        int[] I = MF.RandomPerm(Integers0toNP1);

                        // Find the index in I that is not equal to i and keep the first 3 positions
                        int[] L = MF.RemoveIndex(I,i);
                        int[] r = new int[3];
                        for(int s=0;s<=2;s++)
                            r[s] = L[s];

                        // The three distinct members of the population
                        for(int s=0;s<=4;s++)
                        {
                            Pr1[s] = P[s,r[0]];
                            Pr2[s] = P[s,r[1]];
                            Pr3[s] = P[s,r[2]];
                        }
                        int[] Integers1to5 = { 1,2,3,4,5 };
                        int[] R = MF.RandomPerm(Integers1to5);
                        double[] Pnew = { 0.0,0.0,0.0,0.0,0.0 };

                        // Steps 2 and 3.  Mutation and recombination
                        double Ri;
                        double u;
                        for(int j=0;j<=4;j++)
                        {
                            Ri = R[0];
                            u = U[j,i,k];
                            if((u <= CR) | (j == Ri))
                                Pnew[j] = Pr1[j] + F*(Pr2[j] - Pr3[j]);
                            else
                                Pnew[j] = P0[j];
                        }
                        // Repeat above to code to ensure new members fall within the
                        // range of acceptable parameter values

                        int[] Flag = new int[5] { 0,0,0,0,0 };
                        for(int s=0;s<=4;s++)
                        {
                            if(Pnew[s] <= Lo[s] | Pnew[s] >= Hi[s])
                                Flag[s] = 1;
                        }
                        int Condition = Flag.Sum();
                    if(Condition > 0) goto Step0;

                    // Step 4.  Selection
                    // Calculate the objective function for the ith member and for the candidate
                    double f0 = OF.f(P0,ofset);
                    double fnew = OF.f(Pnew,ofset);

                    // Verify whether the candidate should replace the i-th member
                    // in the population and replace if conditions are satisfied
                    if(fnew < f0)
                        for(int s=0;s<=4;s++) P[s,i] = Pnew[s];
                }
            }

            // Calculate the objective function for each member in the updated population
            double[] fs = new double[NP];
            double[] Pmember = new double[5];
            for(int s=0;s<=NP-1;s++)
            {
                for(int t=0;t<=4;t++)
                {
                    Pmember[t] = P[t,s];
                }
                fs[s] = OF.f(Pmember,ofset);
            }

            // Find the member with the lowest objective function
            double Minf = fs[0];
            int Minpos = 0;
            for(int s=0;s<=NP-1;s++)
            {
                if(fs[s] < Minf)
                {
                    Minf = fs[s];
                    Minpos = s;
                }
            }

            // Return the selected member/parameter
            HParam Pfinal = new HParam();
            Pfinal.kappa = P[0,Minpos];
            Pfinal.theta = P[1,Minpos];
            Pfinal.sigma = P[2,Minpos];
            Pfinal.v0    = P[3,Minpos];
            Pfinal.rho   = P[4,Minpos];

            return Pfinal;
        }
    }
}



