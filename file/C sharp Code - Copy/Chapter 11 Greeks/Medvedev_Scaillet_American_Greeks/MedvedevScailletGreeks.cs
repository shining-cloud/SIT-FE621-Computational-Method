using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Medvedev_Scaillet_American_Greeks
{
    class MSGreeks
    {
        public double MSGreeksFD(HParam param,OpSet opset,int method,double A,double B,int N,double hi,double tol,int MaxIter,int NumTerms,double yinf,string Greek)
        {
            MSExpansionHeston MS = new MSExpansionHeston();
            
            double[] output = new double[6];
            double AmerPut,AmerPutp,AmerPutm;
            double AmerPutpp,AmerPutpm,AmerPutmp,AmerPutmm;

            double S = opset.S;
            double v0 = param.v0;
            double T = opset.T;

            // Define the finite difference increments
            double ds = opset.S  * 0.005;
            double dt = opset.T  * 0.005;
            double dv = param.v0 * 0.005;

            if(Greek == "price")
            {
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPut = output[2];
                return AmerPut;
            }
            else if((Greek == "delta") || (Greek == "gamma"))
            {
                opset.S = S + ds;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutp = output[2];
                opset.S = S - ds;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutm = output[2];
                if(Greek == "gamma")
                {
                    opset.S = S;
                    output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                    AmerPut = output[2];
                    return (AmerPutp - 2.0*AmerPut + AmerPutm)/ds/ds;
                }
                else
                    return (AmerPutp - AmerPutm)/2.0/ds;
            }
            else if(Greek == "theta")
            {
                opset.T = T + dt;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutp = output[2];
                opset.T = T - dt;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutm = output[2];
                return -(AmerPutp - AmerPutm)/2.0/dt;
            }
            else if ((Greek == "vega1") || (Greek == "volga"))
            {
                param.v0 = v0 + dv;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutp = output[2];
                param.v0 = v0 - dv;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutm = output[2];
                double Vega1 = (AmerPutp - AmerPutm)/2.0/dv*2.0*Math.Sqrt(v0);
                if (Greek == "volga")
                {
                    param.v0 = v0;
                    output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                    AmerPut = output[2];
                    double dC2 = (AmerPutp - 2.0*AmerPut + AmerPutm)/dv/dv;
                    return 4.0*Math.Sqrt(v0)*(dC2*Math.Sqrt(v0) + Vega1/4.0/v0);
                }
                else
                    return Vega1;
            }
            else if(Greek == "vanna")
            {
                opset.S  = S + ds;
                param.v0 = v0 + dv;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutpp = output[2];
                param.v0 = v0 - dv;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutpm = output[2];
                opset.S = S - ds;
                param.v0 = v0 + dv;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutmp = output[2];
                param.v0 = v0 - dv;
                output = MS.MSPrice(param,opset,method,A,B,N,hi,tol,MaxIter,NumTerms,yinf);
                AmerPutmm = output[2];
                return (AmerPutpp - AmerPutpm - AmerPutmp + AmerPutmm)/4.0/dv/ds*2.0*Math.Sqrt(v0);
            }
            else
                return 0.0;
        }
    }
}
