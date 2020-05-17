using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Risk_Neutral_Density
{
    class RND
    {
        public fS ExtractRND(double[] K,double[] CallPrice)
        {
            // Calculate the first derivatives of calls w.r.t. strike K. 
            // Use central differences
            int N = K.Length;
            double dK,dC,dC2;
            double[] dCdK = new double[N-2];
            for(int k=1;k<=N-2;k++)
            {
                dK = K[k+1] - K[k-1];
                dC = CallPrice[k+1] - CallPrice[k-1];
                dCdK[k-1] = dC/dK;
            }
            // Calculate the risk neutral density by central finite differences.
            int N2 = dCdK.Length;
            double[] f = new double[N-4];
            for(int k=1;k<=N2-3;k++)
            {
                dK = K[k+1] - K[k-1];
                dC2 = dCdK[k+1] - dCdK[k-1];
                f[k-1] = dC2/dK;
            }
            double[] K2 = new double[N-4];
            for(int k=0;k<=N2-3;k++)
                K2[k] = K[k+2];

            fS output = new fS();
            output.K = K2;
            output.RND = f;
            return output;
        }
    }
}
