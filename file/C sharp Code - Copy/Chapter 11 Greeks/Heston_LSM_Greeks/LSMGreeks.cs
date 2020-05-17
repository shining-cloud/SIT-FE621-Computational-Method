using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Heston_LSM_Greeks
{
    class LSMGreeksAlgo
    {
        public double[] LSMGreeks(OpSet settings,HParam param,int NT,int NS,double[,] Zv,double[,] Zs,string Greek)
        {
            LSM LSM = new LSM();
            Regression R = new Regression();
            MomentMatching MM = new MomentMatching();

            double Spot = settings.S;
            double V0 = param.v0;
            double T = settings.T;
            double r = settings.r;
            double dS = 0.01*Spot;
            double dv = 0.01*V0;
            double dt = 0.01*T;
            double dr = 0.01*r;

            HParam paramP = param;
            HParam paramM = param;
            paramP.v0 = V0+dv;
            paramM.v0 = V0-dv;

            double[] C = new double[2];
            double[] Cp = new double[2];
            double[] Cm = new double[2];
            double[] Cpp = new double[2];
            double[] Cpm = new double[2];
            double[] Cmp = new double[2];
            double[] Cmm = new double[2];
            double Euro=0.0,Amer=0.0;
            double[] output = new double[2]{0.0,0.0};
            if(Greek == "price")
            {
                return LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
            }
            else if((Greek == "delta") | (Greek == "gamma"))
            {
                settings.S = Spot + dS;
                Cp = LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                settings.S = Spot - dS;
                Cm = LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                Euro = (Cp[0] - Cm[0])/2.0/dS;
                Amer = (Cp[1] - Cm[1])/2.0/dS;
                output[0] = Euro;
                output[1] = Amer;
                if(Greek == "delta")
                    return output;
                else
                {
                    C = LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                    Euro = (Cp[0] - 2.0*C[0] + Cm[0])/dS/dS;
                    Amer = (Cp[1] - 2.0*C[1] + Cm[1])/dS/dS;
                    output[0] = Euro;
                    output[1] = Amer;
                    return output;
                }
            }
            else if(Greek == "vega1")
            {
                Cp = LSM.HestonLSM(R.MTrans(MM.MMSim(paramP,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                Cm = LSM.HestonLSM(R.MTrans(MM.MMSim(paramM,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                Euro = (Cp[0] - Cm[0])/2.0/dv*2.0*Math.Sqrt(V0);
                Amer = (Cp[1] - Cm[1])/2.0/dv*2.0*Math.Sqrt(V0);
                output[0] = Euro;
                output[1] = Amer;
                return output;
            }
            else if(Greek == "vanna")
            {
                settings.S = Spot+dt;
                Cpp = LSM.HestonLSM(R.MTrans(MM.MMSim(paramP,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                settings.S = Spot+dt;
                Cpm = LSM.HestonLSM(R.MTrans(MM.MMSim(paramM,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                settings.S = Spot-dt;
                Cmp = LSM.HestonLSM(R.MTrans(MM.MMSim(paramP,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                settings.S = Spot-dt;
                Cmm = LSM.HestonLSM(R.MTrans(MM.MMSim(paramM,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                Euro = (Cpp[0] - Cpm[0] - Cmp[0] + Cmm[0])/2.0/dv/dS*Math.Sqrt(V0);
                Amer = (Cpp[1] - Cpm[1] - Cmp[1] + Cmm[1])/2.0/dv/dS*Math.Sqrt(V0);
                output[0] = Euro;
                output[1] = Amer;
                return output;
            }
            else if(Greek == "theta")
            {
                settings.T = T+dt;
                Cp = LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                settings.T = T-dt;
                Cm = LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                Euro = -(Cp[0] - Cm[0])/2.0/dt;
                Amer = -(Cp[1] - Cm[1])/2.0/dt;
                output[0] = Euro;
                output[1] = Amer;
                return output;
            }
            else if(Greek == "rho")
            {
                settings.r = r+dr;
                Cp = LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                settings.r = r-dr;
                Cm = LSM.HestonLSM(R.MTrans(MM.MMSim(param,settings,NT,NS,Zv,Zs)),settings.K,settings.r,settings.q,settings.T,NT,NS,settings.PutCall);
                Euro = -(Cp[0] - Cm[0])/2.0/dr;
                Amer = -(Cp[1] - Cm[1])/2.0/dr;
                output[0] = Euro;
                output[1] = Amer;
                return output;
            }
            else
                return output;
        }
    }
}

