using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;

namespace Lewis_Vol_of_Vol_Expansion
{
    public partial class MainForm:Form
    {
        public MainForm() { InitializeComponent(); }

        private void GenerateLewisTable_MouseClick(object sender,MouseEventArgs e)
        {
            double S = 100.0;				    // Spot Price
            double T = 0.25;			        // Maturity in Years
            double rf = 0.0;					// Interest Rate
            double q = 0.0;                     // Dividend yield
            string PutCall = "C";               // "P"ut or "Call"
            int trap = 1;                       // 1="Little Trap" characteristic function

            OpSet opset = new OpSet();
            opset.S = S;
            opset.T = T;
            opset.r = rf;
            opset.q = q;
            opset.PutCall = PutCall;
            opset.trap = trap;

            double kappa = 4.0; 			    // Heston Parameter 
            double theta = 0.09/4.0;			// Heston Parameter 
            double sigma = 0.1;				    // Heston Parameter: Volatility of Variance
            double v0 = 0.0225;					// Heston Parameter: Current Variance
            double rho = -0.5;					// Heston Parameter: Correlation
            double lambda = 0.0;				// Heston Parameter 

            HParam param = new HParam();
            param.kappa = kappa;
            param.theta = theta;
            param.sigma = sigma;
            param.v0 = v0;
            param.rho = rho;
            param.lambda = lambda;

            int NK = 7;
            double[] K = new double[7]          // Strikes
            { 70.0,80.0,90.0,100.0,110.0,120.0,130.0 };

            // Time average of the deterministic variance and volatility
            double vol = param.theta + (param.v0 - param.theta)*(1.0-Math.Exp(-param.kappa*T))/(param.kappa*T);
            double IV = Math.Sqrt(vol);

            // Classes
            BlackScholes BS = new BlackScholes();
            NewtonCotes NC = new NewtonCotes();
            Lewis L = new Lewis();
            BisectionAlgo BA = new BisectionAlgo();

            // Table 3.3.1 on Page 81 of Lewis (2001)
            double[] SeriesIPrice = new double[NK];         // Series I price and implied vol
            double[] IVI = new double[NK];
            double[] SeriesIIPrice = new double[NK];        // Series II price and implied vol
            double[] IVII = new double[NK];
            double[] HPrice = new double[NK];                // Exact Heston price and implied vol
            double[] IVe = new double[NK];
            double[] BSPrice = new double[NK];              // Black Scholes price

            // Bisection algorithm settings
            double lo = 0.01;
            double hi = 2.0;
            double Tol = 1.0e-10;
            int MaxIter = 10000;

            // Newton Cotes settings
            int method = 3;
            double a = 1.0e-10;
            double b = 100.0;
            int N = 10000;

            double[] SeriesII = new double[2];
            for(int k=0;k<=NK-1;k++)
            {
                // Generate the prices and implied vols
                opset.K = K[k];
                SeriesIPrice[k] = L.SeriesICall(S,K[k],rf,q,T,v0,rho,theta,kappa,sigma);
                IVI[k] = BA.BisecBSIV(PutCall,S,K[k],rf,q,T,lo,hi,SeriesIPrice[k],Tol,MaxIter);
                SeriesII = L.SeriesIICall(S,K[k],rf,q,T,v0,rho,theta,kappa,sigma);
                SeriesIIPrice[k] = SeriesII[0];
                IVII[k] = SeriesII[1];
                HPrice[k] = NC.HestonPriceNewtonCotes(param,opset,method,a,b,N);
                IVe[k] = BA.BisecBSIV(PutCall,S,K[k],rf,q,T,lo,hi,HPrice[k],Tol,MaxIter);
                BSPrice[k] = BS.BSC(S,K[k],rf,q,vol,T);

                // Write to the List Boxes
                Strike.Items.Add(K[k]);
                ExactPrice.Items.Add(Math.Round(HPrice[k],6));
                IVExact.Items.Add(Math.Round(100*IVe[k],2));
                SeriesICall.Items.Add(Math.Round(SeriesIPrice[k],6));
                IV1.Items.Add(Math.Round(100*IVI[k],2));
                SeriesIICall.Items.Add(Math.Round(SeriesIIPrice[k],6));
                IV2.Items.Add(Math.Round(100*IVII[k],2));
                BlackScholesPrice.Items.Add(Math.Round(BSPrice[k],6));
                IVBS.Items.Add(Math.Round(100*IV,2));
            }
        }

        private void label2_Click(object sender,EventArgs e)
        {

        }
    }
}
