using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Gauthier_Starting_Values
{
    class ObjectiveFunction
    {
        // Returns the objective function for calculating the Gauthier and Rivaille coefficients under optimization
        public double f(double[] param,double[] Coeff1,double[] Coeff2,double Put1,double Put2)
        {
            double sigma = param[0];
            double rho   = param[1];

            double A1 = Coeff1[0];
            double B1 = Coeff1[1];
            double C1 = Coeff1[2];
            double D1 = Coeff1[3];

            double A2 = Coeff2[0];
            double B2 = Coeff2[1];
            double C2 = Coeff2[2];
            double D2 = Coeff2[3];

            double ObjFun = Math.Pow(A1 + B1*sigma*sigma + C1*rho*sigma + D1*rho*rho*sigma*sigma - Put1,2) + 
                            Math.Pow(A2 + B2*sigma*sigma + C2*rho*sigma + D2*rho*rho*sigma*sigma - Put2,2);
            return ObjFun;
        }
    }
}
 
