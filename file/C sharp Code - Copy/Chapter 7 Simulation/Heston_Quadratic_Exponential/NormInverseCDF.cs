using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using System.IO;

namespace Heston_Quadratic_Exponential
{
    class InverseNorm
    {
        public double normICDF(double p)
        {
            // Calculates the inverse standard normal CDF

            // Coefficients for p close to 0.5
            double a0 = 3.3871328727963666080e0;
            double a1 = 1.3314166789178437745e+2;
            double a2 = 1.9715909503065514427e+3;
            double a3 = 1.3731693765509461125e+4;
            double a4 = 4.5921953931549871457e+4;
            double a5 = 6.7265770927008700853e+4;
            double a6 = 3.3430575583588128105e+4;
            double a7 = 2.5090809287301226727e+3;
            double b1 = 4.2313330701600911252e+1;
            double b2 = 6.8718700749205790830e+2;
            double b3 = 5.3941960214247511077e+3;
            double b4 = 2.1213794301586595867e+4;
            double b5 = 3.9307895800092710610e+4;
            double b6 = 2.8729085735721942674e+4;
            double b7 = 5.2264952788528545610e+3;

            // Coefficients for p not close to 0, 0.5, or 1
            double c0 = 1.42343711074968357734e0;
            double c1 = 4.63033784615654529590e0;
            double c2 = 5.76949722146069140550e0;
            double c3 = 3.64784832476320460504e0;
            double c4 = 1.27045825245236838258e0;
            double c5 = 2.41780725177450611770e-1;
            double c6 = 2.27238449892691845833e-2;
            double c7 = 7.74545014278341407640e-4;
            double d1 = 2.05319162663775882187e0;
            double d2 = 1.67638483018380384940e0;
            double d3 = 6.89767334985100004550e-1;
            double d4 = 1.48103976427480074590e-1;
            double d5 = 1.51986665636164571966e-2;
            double d6 = 5.47593808499534494600e-4;
            double d7 = 1.05075007164441684324e-9;

            // Coefficients for p near 0 or 1
            double e0 = 6.65790464350110377720e0;
            double e1 = 5.46378491116411436990e0;
            double e2 = 1.78482653991729133580e0;
            double e3 = 2.96560571828504891230e-1;
            double e4 = 2.65321895265761230930e-2;
            double e5 = 1.24266094738807843860e-3;
            double e6 = 2.71155556874348757815e-5;
            double e7 = 2.01033439929228813265e-7;
            double f1 = 5.99832206555887937690e-1;
            double f2 = 1.36929880922735805310e-1;
            double f3 = 1.48753612908506148525e-2;
            double f4 = 7.86869131145613259100e-4;
            double f5 = 1.84631831751005468180e-5;
            double f6 = 1.42151175831644588870e-7;
            double f7 = 2.04426310338993978564e-15;

            double q = p - 0.5;
            double y=0.0;
            double r=0.0;
            if(Math.Abs(q) < 0.425)
            {
                // For p close to 0.5
                r = 0.180625 - q*q;
                y = q*(((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / 
                    (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + 1);
            }
            else
            {
                if(q < 0.0)
                    r = p;
                else
                    r = 1.0 - p;
                if(r <= 0.0)
                    y = 0.0;
                r = Math.Sqrt(-Math.Log(r));
                if(r <= 5.0)
                {
                    // For p not close to 0, 0.5, or 1
                    r = r - 1.6;
                    y = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) /
                        (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + 1);
                }
                else
                {
                    // For p near 0 or 1
                    r = r - 5.0;
                    y = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) /
                        (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + 1);
                }
                if(q < 0.0)
                    y = -y;
            }
            return y;
        }
    }
}


