using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Beliaeva_Nawalkha_Tree
{
    class Probabilities
    {
        public float[] probV(float X,float X0,float dt,float kappa,float theta,float sigma)
        {
            // Obtain the V-tree probabilities for the Belieava-Nawalkha tree
            // INPUTS
            //   X = single transformed volatility
            //   X0 = time zero transformed volatility
            //   dt = time increment
            //   k  = ceiling(Math.Sqrt(V(t)/V(0))
            //   rho = Heston CIR parameter
            //   sigma = Heston CIR parameter
            //   kappa = Heston CIR parameter
            // OUTPUTS
            //  Matrix of probabilities pu,pm,pd

            // Equations (24) through (26)
            float be = Convert.ToSingle(X0/Math.Sqrt(dt)/Math.Floor(X0/Math.Sqrt(1.5*dt)));
            float bc = Convert.ToSingle(X0/Math.Sqrt(dt)/Math.Floor(X0/Math.Sqrt(1.5*dt)+1.0));
            float b;
            if(Math.Abs(bc-Math.Sqrt(1.5)) < Math.Abs(be-Math.Sqrt(1.5)))
                b = bc;
            else
                b = be;

            // Equations (22) and (30)
            float muX = Convert.ToSingle(1.0/X*(0.5*kappa*(4.0*theta/sigma/sigma-X*X)-0.5));
            float J =  Convert.ToSingle(Math.Floor(muX*Math.Sqrt(dt)/b + 1/b/b));
            float pu = 0.0f;
            float pm = 0.0f;
            float pd = 0.0f;

            if(X > 0)
            {
                // Probabilities where X > 0 (Equation 28)
                pu = Convert.ToSingle(0.5/b/b - 0.5*J + 0.5/b*muX*Math.Sqrt(dt));
                pm = 1.0f - 1.0f/b/b;
                pd = Convert.ToSingle(0.5/b/b + 0.5*J - 0.5/b*muX*Math.Sqrt(dt));
            }
            else
            {
                // Probabilities where X = 0 (Equation 33)
                float Xu = Convert.ToSingle(X + b*(J+1.0)*Math.Sqrt(dt));
                float Vu = Xu*Xu*sigma*sigma/4.0f;
                pu = kappa*theta*dt/Vu;
                pm = 0.0f;
                pd = 1.0f - pu;
            }

            // Output the results
            float[] output = new float[3];
            output[0] = pu;
            output[1] = pm;
            output[2] = pd;
            return output;
        }
        public float[] probY(float Vt,float V0,float Yt,float dt,float rho,float sigma,float kappa)
        {
            // Obtain the Y-tree probabilities for the Belieava-Nawalkha tree
            // INPUTS
            //   Vt = current value of V(t)
            //   V0 = value of V(0) (Heston parameter)
            //   dt = time increment
            //   k  = ceiling(Math.Sqrt(V(t)/V(0))
            //   rho = Heston CIR parameter
            //   sigma = Heston CIR parameter
            //   kappa = Heston CIR parameter
            // OUTPUTS
            //   Single probabilities pu,pm,pd

            // Equation 11
            float k;
            if(Vt > 0.0)
                k = Convert.ToSingle(Math.Ceiling(Math.Sqrt(Vt/V0)));
            else
                k = 1.0f;

            // Equations 9 and 6
            float muY = Convert.ToSingle((rho/sigma*kappa - 0.5)*Vt);
            float sigmayt = Convert.ToSingle(Math.Sqrt(1-rho*rho)*Math.Sqrt(Vt));
            float sigmay0 = Convert.ToSingle(Math.Sqrt(1-rho*rho)*Math.Sqrt(V0));

            // Calculate the Yu, Ym, Yd
            float I = Convert.ToSingle((Math.Round(muY/k/sigmay0*Math.Sqrt(dt))));
            float Yu = Convert.ToSingle(Yt + (I+1.0)*k*sigmay0*Math.Sqrt(dt));
            float Ym = Convert.ToSingle(Yt + (I+0.0)*k*sigmay0*Math.Sqrt(dt));
            float Yd = Convert.ToSingle(Yt + (I-1.0)*k*sigmay0*Math.Sqrt(dt));

            // Equation 17
            float eu = Yu - Yt - muY*dt;
            float em = Ym - Yt - muY*dt;
            float ed = Yd - Yt - muY*dt;

            // Equation 16
            float pu = Convert.ToSingle(0.5*(sigmayt*sigmayt*dt + em*ed) / (k*k*sigmay0*sigmay0) / dt);
            float pm = Convert.ToSingle(   -(sigmayt*sigmayt*dt + eu*ed) / (k*k*sigmay0*sigmay0) / dt);
            float pd = Convert.ToSingle(0.5*(sigmayt*sigmayt*dt + eu*em) / (k*k*sigmay0*sigmay0) / dt);

            // Output the results
            float[] output = new float[3];
            output[0] = pu;
            output[1] = pm;
            output[2] = pd;
            return output;
        }
    }
}
