using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace Beliaeva_Nawalkha_Tree
{
    class MainTree
    {
        float[,] X;
        float[,] V;

        public BivariateStruct BuildBivariateTree(float S0,string PutCall,float Strike,float T,float rf,int NT,float kappa,float theta,float sigma,float V0,float rho,float threshold)
        {
            float sigmay0 = Convert.ToSingle(Math.Sqrt(1-rho*rho)*Math.Sqrt(V0));
            float dt = T/NT;

            // Generate the volatility trees
            VolTree VT = new VolTree();
            VolStruct VolTree = VT.BuildVolTree(kappa,theta,sigma,V0,dt,NT,threshold);
            X = VolTree.X;
            V = VolTree.V;

            // Erase the VolTree structured object since it's no longer needed
            VolTree.Dispose();

            int RBound = VolTree.RBound;
            int M = VolTree.M;

//            if(RBound != 2*NT-1)
//                Console.WriteLine("Vt matrix not triangular");

            // Find the column where the tree changes from triangular to rectangular
            int ColChange = RBound - M + 1;

            // The values of k(t) from Equation (11) and I(t) from (14)
            int[,] k = new int[2*NT-1,NT];
            for(int t=0;t<=NT-1;t++)
                for(int n=M-t;n<=M+t;n++)
                {
                    if(V[n,t] > 0)
                        k[n,t] = Convert.ToInt16(Math.Ceiling(Math.Sqrt(V[n,t]/V0)));
                    else
                        k[n,t] = 1;
                }

            // Break up k[n,t] into columns and store the maximum of each column
            int[] maxK = new int[NT];
            int[] kcol = new int[2*NT-1];

            for(int t=0;t<=NT-1;t++)
            {
                for(int n=0;n<=2*NT-2;n++)
                    kcol[n] = k[n,t];
                maxK[t] = kcol.Max();
            }

            int[] numY = new int[NT];
            int[] numV = new int[NT];
            int[] numRows = new int[NT];
            numY[0] = 1;
            numY[1] = 3;
            numV[0] = 1;
            numV[1] = 3;
            numRows[0] = 1;
            numRows[1] = 9;
            for(int t=2;t<=NT-1;t++)
            {
                numY[t] = numY[t-1] + 2*maxK[t];
                numV[t] = 2*(t+1) - 1;
                numRows[t] = numV[t]*numY[t];
            }
            // Number of rows required for the main Yt(n,t) matrix
            int NR = numRows.Max();
            Console.WriteLine("Yt matrix for log stock price has {0:0,0} rows and {1:0,0} columns = {2:0,0} elements",NR,NT,NR*NT);

            // Find the branch indices
            int[,] Branch = new int[numRows[NT-2],9*(NT-1)];
            for(int j=0;j<=8;j++)
                Branch[0,j] = j;

            for(int t=1;t<=NT-2;t++)
            {
                int nY = numY[t];
                int First = maxK[t];
                int nK = 2*t+1;
                int[] K = new int[nK];
                for(int j=0;j<=nK-1;j++)
                    K[j] = k[M-t+j,t];
                for(int n=0;n<=numRows[t]-1;n++)
                {
                    int a = Convert.ToInt32(Math.Ceiling(Convert.ToDouble(n/nY)));
                    int b = n % nY;

                    // Find the middle-to-middle branches
                    Branch[n,9*t+1] = First + a*numY[t+1] + b;
                    Branch[n,9*t+4] = Branch[n,9*t+1] +   numY[t+1];
                    Branch[n,9*t+7] = Branch[n,9*t+1] + 2*numY[t+1];

                    // Find the rest of the branches
                    Branch[n,9*t+0] = Branch[n,9*t+1] - K[a];
                    Branch[n,9*t+2] = Branch[n,9*t+1] + K[a];
                    Branch[n,9*t+3] = Branch[n,9*t+4] - K[a];
                    Branch[n,9*t+5] = Branch[n,9*t+4] + K[a];
                    Branch[n,9*t+6] = Branch[n,9*t+7] - K[a];
                    Branch[n,9*t+8] = Branch[n,9*t+7] + K[a];
                }
            }

            // Adjust the last branch upward
            for(int t=ColChange;t<=NT-2;t++)
                for(int j=6;j<=8;j++)
                    Branch[numRows[t]-1,9*t+j] = Branch[numRows[t]-2,9*t+j];

            // Find the values for the indices and for the probabilities
            float Vt = new float();
            float Xt = new float();
            int Kt = new int();

            // Log stock tree (Yt), stock price tree (St), and probabilities (Prob)
            float[,] Yt = new float[NR,NT];
            float Y0 = Convert.ToSingle(Math.Log(S0) - rho*V0/sigma);
            Yt[0,0] = Y0;
            float[,] Prob = new float[numRows[NT-2],9*(NT-1)];
            float X0 = X[M,0];

            Probabilities PR = new Probabilities(); 
            float ht,muy;
            int I;
            int[] NewBranch = new int[9];
            for(int t=0;t<=NT-2;t++)
            {
                int n = -1;
                int[] J = new int[2*t+1];
                for(int j=0;j<=2*t;j++) J[j] = -t + j;
                for(int j=0;j<=numV[t]-1;j++)
                {
                    for(int r=0;r<=numY[t]-1;r++)
                    {
                        n += 1;
                        Vt = V[M+J[j],t];   // Volatility
                        Xt = X[M+J[j],t];   // Transformed volatility
                        Kt = k[M+J[j],t];   // Node jump k(t) 
                        for(int s=0;s<=8;s++)
                            NewBranch[s] = Branch[n,9*t+s];
                        muy = Convert.ToSingle((rho*kappa/sigma - 0.5)*Vt);
                        I = Convert.ToInt32(Math.Round(muy/Kt/sigmay0*Math.Sqrt(dt)));
                        if(Yt[n,t] > 0)
                        {
                            for(int s=0;s<=8;s++)
                            {
                                if(s==1 || s==4 || s==7)
                                    // Middle node
                                    Yt[NewBranch[s],t+1] = Yt[n,t] + (I+0)*Kt*sigmay0*Convert.ToSingle(Math.Sqrt(dt));
                                else if(s==0 || s==3 || s==6)
                                    // Up node
                                    Yt[NewBranch[s],t+1] = Yt[n,t] + (I+1)*Kt*sigmay0*Convert.ToSingle(Math.Sqrt(dt));
                                else if(s==2 || s==5 || s==8)
                                    // Down node
                                    Yt[NewBranch[s],t+1] = Yt[n,t] + (I-1)*Kt*sigmay0*Convert.ToSingle(Math.Sqrt(dt));
                            }
                        }
                        if(Yt[n,t] > 0.0f)  // Construct the joint probabilities
                        {
                            float[] pv = PR.probV(Xt,X0,dt,kappa,theta,sigma);
                            float pvu = pv[0];
                            float pvm = pv[1];
                            float pvd = pv[2];
                            float[] py = PR.probY(Vt,V0,Yt[n,t],dt,rho,sigma,kappa);
                            float pyu = py[0];
                            float pym = py[1];
                            float pyd = py[2];
                            float[] prob = new float[9];
                            prob[0] = pvu*pyu;
                            prob[1] = pvu*pym;
                            prob[2] = pvu*pyd;
                            prob[3] = pvm*pyu;
                            prob[4] = pvm*pym;
                            prob[5] = pvm*pyd;
                            prob[6] = pvd*pyu;
                            prob[7] = pvd*pym;
                            prob[8] = pvd*pyd;
                            for(int s=0;s<=8;s++)
                                Prob[n,9*t+s] = prob[s];
                        }
                    }
                }
            }

            // Find the American and European option prices
            float[,] Euro = new float[NR,NT];
            float[,] Amer = new float[NR,NT];

            // Stock price at maturity
            float[] ST = new float[numRows[NT-1]];
            ht = Convert.ToSingle((rf - rho*kappa*theta/sigma)*Convert.ToDouble(NT-1)*dt);
            int[] JJ = new int[2*(NT-1)+1];
            for(int j=0;j<=2*(NT-1);j++)
                JJ[j] = -(NT-1) + j;
            int m = -1;
            for(int j=0;j<=numV[NT-1]-1;j++)
            {
                for(int r=0;r<=numY[NT-1]-1;r++)
                {
                    m += 1;
                    Vt = V[M+JJ[j],NT-1];
                    if(Yt[m,NT-1] > 0)
                        ST[m] = Convert.ToSingle(Math.Exp(Yt[m,NT-1] + rho/sigma*Vt + ht));
                }
            }

            // Payoff at maturity
            for(int n=0;n<=NR-1;n++)
            {
                if(PutCall == "C")
                {
                    Euro[n,NT-1] = Convert.ToSingle(Math.Max(ST[n] - Strike,0.0));
                    Amer[n,NT-1] = Convert.ToSingle(Math.Max(ST[n] - Strike,0.0));
                }
                else if(PutCall == "P")
                {
                    Euro[n,NT-1] = Convert.ToSingle(Math.Max(Strike - ST[n],0.0));
                    Amer[n,NT-1] = Convert.ToSingle(Math.Max(Strike - ST[n],0.0));
                }
            }

            float[] P = new float[9];
            int[] B = new int[9];
            float St = new float();

            for(int t=NT-2;t>=0;t--)
            {
                int n = -1;
                ht = Convert.ToSingle((rf - rho*kappa*theta/sigma)*Convert.ToDouble(t)*dt);
                int[] J = new int[2*t+1];
                for(int j=0;j<=2*t;j++)
                    J[j] = -t + j;
                for(int j=0;j<=numV[t]-1;j++)
                {
                    for(int r=0;r<=numY[t]-1;r++)
                    {
                        n += 1;
                        Vt = V[M+J[j],t];
                        if(Yt[n,t] > 0)
                        {
                            St = Convert.ToSingle(Math.Exp(Yt[n,t] + rho/sigma*Vt + ht));  // Stock price
                            for(int s=0;s<=8;s++)
                            {
                                P[s] = Prob[n,9*t+s];
                                B[s] = Branch[n,9*t+s];
                            }
                            Euro[n,t] = P[0]*Euro[B[0],t+1] + P[1]*Euro[B[1],t+1] + P[2]*Euro[B[2],t+1]   
                                  + P[3]*Euro[B[3],t+1] + P[4]*Euro[B[4],t+1] + P[5]*Euro[B[5],t+1]   
                                  + P[6]*Euro[B[6],t+1] + P[7]*Euro[B[7],t+1] + P[8]*Euro[B[8],t+1];
                            Euro[n,t] = Convert.ToSingle(Math.Exp(-rf*dt)*Euro[n,t]);
                            Amer[n,t] = P[0]*Amer[B[0],t+1] + P[1]*Amer[B[1],t+1] + P[2]*Amer[B[2],t+1]   
                                  + P[3]*Amer[B[3],t+1] + P[4]*Amer[B[4],t+1] + P[5]*Amer[B[5],t+1]   
                                  + P[6]*Amer[B[6],t+1] + P[7]*Amer[B[7],t+1] + P[8]*Amer[B[8],t+1];
                            Amer[n,t] = Convert.ToSingle(Math.Exp(-rf*dt)*Amer[n,t]);
                            if(PutCall == "C")
                                Amer[n,t] = Math.Max(St - Strike,Amer[n,t]);
                            else if(PutCall == "P")
                                Amer[n,t] = Math.Max(Strike - St,Amer[n,t]);
                        }
                    }
                }
            }
            float EuroPrice = Euro[0,0];  // European option
            float AmerPrice = Amer[0,0];  // American option

            // Output the results
            BivariateStruct output = new BivariateStruct();
            output.Euro = Euro;
            output.Amer = Amer;
            output.Yt = Yt;
            output.V = V;
            output.X = X;
            //            output.St = St;
            output.Prob = Prob;
            output.EuroPrice = EuroPrice;
            output.AmerPrice = AmerPrice;
            return output;
        }
    }
}

