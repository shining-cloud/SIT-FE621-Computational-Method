using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace ADI_Method
{
    class ADIMethod
    {
        public double ADIPrice(string scheme,double thet,HParam param,double S0,double V0,double K,double r,double q,double[] S,double[] V,double[] T,string GridType)
        {
            // Alternating Direction Implicit (ADI) scheme for the Heston model.
            // INPUTS
            //   scheme = 'DO' Douglas, 'CS' Craig-Sneyd, 
            //            'MSC' Modified CS, 'HV' Hundsdorfer-Verwer
            //   thet  = weighing parameter 0 = Explicit Scheme
            //           1 = Implicit Scheme, 0.5 = Crank-Nicolson
            //  param = Heston parameters
            //  S0 = Spot price on which to price
            //  V0 = Volatility on which to price
            //  K = Strike price
            //  r = Risk free rate
            //  q = Dividend yield
            //  S = Stock price grid - uniform
            //  V = Volatility grid - uniform
            //  T = Maturity grid - uniform

            Interpolation IP = new Interpolation();
            MatrixOps MO = new MatrixOps();

            // Heston parameters
            double kappa = param.kappa;
            double theta = param.theta;
            double sigma = param.sigma;
            double v0    = param.v0;
            double rho   = param.rho;
            double lambda = param.lambda;

            // Length of grids;
            int NS = S.Length;
            int NV = V.Length;
            int NT = T.Length;
            double dt = T[1]-T[0];

            // Build the matrices for the derivatives and other matrices
            BuildDerivativesNU BN = new BuildDerivativesNU();
            BuildDerivativesU BU = new BuildDerivativesU();
            LMatrices LMat;
            if(GridType == "Uniform")
                LMat = BU.BuildDerivatives(S,V,T);
            else
                LMat = BN.BuildDerivativesNonUniform(S,V,T);

            double[,] derS  = LMat.derS;
            double[,] derSS = LMat.derSS;
            double[,] derV1 = LMat.derV1;
            double[,] derV2 = LMat.derV2;
            double[,] derVV = LMat.derVV;
            double[,] derSV = LMat.derSV;
            double[,] R     = LMat.R;

            // Decompose the derivatives matrices and create the identity matrix
            int N = NS*NV;
            double[,] A0 = new double[N,N];
            double[,] A1 = new double[N,N];
            double[,] A2 = new double[N,N];
            double[,] I  = new double[N,N];
            for(int i=0;i<=N-1;i++)
            {
                I[i,i] = 1.0;
                for(int j=0;j<=N-1;j++)
                {
                    A0[i,j] = rho*sigma*derSV[i,j];
                    A1[i,j] = (r-q)*derS[i,j] + 0.5*derSS[i,j] - 0.5*r*R[i,j];
                    A2[i,j] = kappa*theta*derV1[i,j] - kappa*derV2[i,j] + 0.5*sigma*sigma*derVV[i,j] - 0.5*r*R[i,j];
                }
            }

            // Initialize the u vector, create identity matrix
            // u plays the role of U(t-1)
            double[] U = new double[N];
            double[] u = new double[N];

            // U(0) vector - value of U(T) at maturity
            double[] Si = new double[N];
            int k = 0;
            for(int v=0;v<=NV-1;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    Si[k] = S[s];
                    U[k] = Math.Max(Si[k] - K,0.0);
                    k += 1;
                }

            // Matrices for the ADI method
            double[] Y0 = new double[N];
            double[] Y1 = new double[N];
            double[] Y2 = new double[N];
            double[] Y0_ = new double[N];
            double[] Y1_ = new double[N];
            double[] Y2_ = new double[N];
            double[] Y0h = new double[N];

            //// Loop through the time increment
            for(int t=1;t<=NT-1;t++)
            {
                u  = U;

                // Vectors common to all ADI schemes
                double[,] SumA = MO.MAdd(MO.MAdd(A0,A1),A2);
                double[,] M1   = MO.MMultS(SumA,dt);
                double[,] IM1  = MO.MAdd(I,M1);
                Y0  = MO.MVMult(IM1,u);

                double[,] dtA1 = MO.MMultS(A1,dt*thet);
                double[,] IminusA1 = MO.MSub(I,dtA1);
                double[,] A1dt = MO.MMultS(A1,dt*thet);
                double[] A1u = MO.MVMult(A1dt,u);
                double[] Y0minusA1 = MO.VSub(Y0,A1u);
                Y1 = MO.MVMult(MO.MInvLU(IminusA1),Y0minusA1);

                double[,] dtA2 = MO.MMultS(A2,dt*thet);
                double[,] IminusA2 = MO.MSub(I,dtA2);
                double[,] A2dt = MO.MMultS(A2,dt*thet);
                double[] A2u = MO.MVMult(A2dt,u);
                double[] Y1minusA2 = MO.VSub(Y1,A2u);
                Y2 = MO.MVMult(MO.MInvLU(IminusA2),Y1minusA2);

                if(scheme == "DO")
                {
                    // Douglas ADI scheme
                    U = Y2;
                }
                else if(scheme == "CS")
                {
                    // Craig-Sneyd ADI scheme
                    double[] A0Y2 = MO.MVMult(A0,Y2);
                    double[] A0u  = MO.MVMult(A0,u);
                    double[] A0A0 = MO.VSub(A0Y2,A0u);
                    double[] dtA  = MO.VMultS(A0A0,dt*0.5);
                    Y0_ = MO.VAdd(Y0,dtA);

                    double[] Y0_minusA1 = MO.VSub(Y0_,A1u);
                    Y1_ = MO.MVMult(MO.MInvLU(IminusA1),Y0_minusA1);

                    double[] Y1_minusA2 = MO.VSub(Y1_,A2u);
                    Y2_ = MO.MVMult(MO.MInvLU(IminusA2),Y1_minusA2);

                    U = Y2_;
                }
                else if(scheme == "MCS")
                {
                    // Modified Craig-Sneyd ADI scheme
                    double[] A0Y2 = MO.MVMult(A0,Y2);
                    double[] A0u  = MO.MVMult(A0,u);
                    double[] A0A0 = MO.VSub(A0Y2,A0u);
                    double[] dtA  = MO.VMultS(A0A0,dt*thet);
                    Y0h = MO.VAdd(Y0,dtA);

                    double[] SumAY2 = MO.MVMult(SumA,Y2);
                    double[] SumAu  = MO.MVMult(SumA,u);
                    double[] ThetDt = MO.VMultS(MO.VSub(SumAY2,SumAu),(0.5-thet)*dt);
                    Y0_ = MO.VAdd(Y0h,ThetDt);

                    double[] Y0_minusA1 = MO.VSub(Y0_,A1u);
                    Y1_ = MO.MVMult(MO.MInvLU(IminusA1),Y0_minusA1);

                    double[] Y1_minusA2 = MO.VSub(Y1_,A2u);
                    Y2_ = MO.MVMult(MO.MInvLU(IminusA2),Y1_minusA2);

                    U = Y2_;
                }
                else if(scheme == "HV")
                {
                    // Hundsdorfer-Verver ADI scheme
                    double[] A0Y2 = MO.MVMult(SumA,Y2);
                    double[] A0u  = MO.MVMult(SumA,u);
                    double[] A0A0 = MO.VSub(A0Y2,A0u);
                    double[] dtA  = MO.VMultS(A0A0,dt*0.5);
                    Y0_ = MO.VAdd(Y0,dtA);

                    double[] A1Y2 = MO.VMultS(MO.MVMult(A1,Y2),thet*dt);
                    double[] Y0_minusA1Y2 = MO.VSub(Y0_,A1Y2);
                    Y1_ = MO.MVMult(MO.MInvLU(IminusA1),Y0_minusA1Y2);

                    double[] A2Y2 = MO.VMultS(MO.MVMult(A2,Y2),thet*dt);
                    double[] Y1_minusA2Y2 = MO.VSub(Y1_,A2Y2);
                    Y2_ = MO.MVMult(MO.MInvLU(IminusA2),Y1_minusA2Y2);

                    U = Y2_;
                }
            }
            // Restack the U vector to output a matrix
            double[,] UU = new double[NS,NV];
            k = 0;
            for(int v=0;v<=NV-1;v++)
                for(int s=0;s<=NS-1;s++)
                {
                    UU[s,v] = U[k];
                    k += 1;
                }

            // Interpolate to get the price at S0 and v0
            return IP.interp2(V,S,UU,V0,S0);
        }
    }
}
