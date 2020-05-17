using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using Excel = Microsoft.Office.Interop.Excel;

namespace Implied_Volatility
{
    public partial class Form1:Form
    {
        double[,] ModelPrice;
        double[,] ModelIV;
        double[,] MktPrice;
        double[,] MktIV;
        double[] K;
        double[] T;
        int NK,NT;

        public Form1() { InitializeComponent(); }

        public void GeneratePrices_MouseClick_1(object sender,MouseEventArgs e)
        {
            // Gauss Laguerre 32 abscissas and weights
            double[] x = new Double[32];
            double[] w = new Double[32];
            using(TextReader reader = File.OpenText("../../GaussLaguerre32.txt"))
                for(int k=0;k<=31;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    x[k] = double.Parse(bits[0]);
                    w[k] = double.Parse(bits[1]);
                }

            // Size of strikes (NK) and maturities (NT)
            NK = 27;
            NT = 4;

            // Market prices of SPX puts
            string PutCall = "P";
            MktPrice = new Double[NK,NT];
            using(TextReader reader = File.OpenText("../../SPX_MktPrice.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    MktPrice[k,0] = double.Parse(bits[0]);
                    MktPrice[k,1] = double.Parse(bits[1]);
                    MktPrice[k,2] = double.Parse(bits[2]);
                    MktPrice[k,3] = double.Parse(bits[3]);
                }

            // SPX strikes
            K = new Double[NK];
            using(TextReader reader = File.OpenText("../../SPX_K.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(' ');
                    K[k] = double.Parse(bits[0]);
                }

            // SPX maturities
            T = new Double[4] { 0.019178082191781,0.041095890410959,0.117808219178082,0.194520547945205 };
            for(int t=0;t<=NT-1;t++)
                T[t] = Math.Round(T[t],3);

            // Settings for the bisection algorithm
            double a = 0.001;
            double b = 2.5;
            double Tol = 1e-6;
            int MaxIter = 5000;

            // SPX spot price, risk free rate, and dividend yield
            double S = 1164.97;
            double r = 0.0;
            double q = 0.0;
            double trap = 1;

            // Calculate the SPX implied volatilities and write to a text file
            BisectionAlgorithm B = new BisectionAlgorithm();
            MktIV = new Double[NK,NT];
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<=NT-1;t++)
                    MktIV[k,t] = B.BisecBSIV(PutCall,S,K[k],r,q,T[t],a,b,MktPrice[k,t],Tol,MaxIter);

            // Write the SPX IV to a text file
            using(var writer = new StreamWriter("../../SPX_Implied_Volatilities.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    for(int t=0;t<=NT-1;t++)
                    {
                        writer.Write(MktIV[k,t]);
                        writer.Write(' ');
                    }
                    writer.WriteLine();
                }

            // Input the pre-calculated Heston parameters
            HParam param = new HParam();
            param.v0     =  0.12581983608881;
            param.theta  =  0.00000000008998;
            param.kappa  =  0.02271544261755;
            param.sigma  =  0.77411910754829;
            param.rho    = -0.95555303183327;
            param.lambda = 0.0;

            OpSet settings = new OpSet();
            settings.S = S;
            settings.r = r;
            settings.q = q;
            settings.trap = 1;

            // Fit the Heston Model Prices and Implied Vols
            HestonPrice HP = new HestonPrice();
            ModelPrice = new Double[NK,NT];
            ModelIV = new Double[NK,NT];
            for(int k=0;k<=NK-1;k++)
            {
                settings.K = K[k];
                for(int t=0;t<=NT-1;t++)
                {
                    settings.T = T[t];
                    ModelPrice[k,t] = HP.HestonPriceGaussLaguerre(param,settings,x,w);
                    ModelPrice[k,t] = Math.Round(ModelPrice[k,t],2);
                    ModelIV[k,t] = B.BisecBSIV(PutCall,S,K[k],r,q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                }
            }

            // Write the Heston model price to a text file
            using(var writer = new StreamWriter("../../SPX_Model_Price.txt"))
                for(int k=0;k<=NK-1;k++)
                {
                    for(int t=0;t<=NT-1;t++)
                    {
                        writer.Write(ModelPrice[k,t]);
                        writer.Write(' ');
                    }
                    writer.WriteLine();
                }

            // Term structure of ATM variance.  Spot is 1165.97 so take K=1165 as ATM
            int atm = 14;
            double[] MktVol = new Double[NT];
            for(int t=0;t<=NT-1;t++)
                MktVol[t] = MktIV[atm,t] * Math.Sqrt(T[t]);

            // Prime parameters
            double kappa = param.kappa;
            double rho   = param.rho;
            double sigma = param.sigma;
            double theta = param.theta;
            double v0    = param.v0;
            double kappa_ = kappa - rho*sigma/2;
            double theta_ = theta*kappa/kappa_;

            // ATM variance and volatility, Equation (3.18) Gatheral
            double[] ModelVar = new Double[200];
            double[] ModelVol = new Double[200];
            double Time = 0.0001;
            double increment = 0.001;
            for(int t=0;t<=199;t++)
            {
                if(t>0) Time = Time + increment;
                ModelVar[t] = (v0 - theta_) * (1 - Math.Exp(-kappa_*Time)) / (kappa_*Time) + theta_;
                ModelVol[t] = Math.Sqrt(ModelVar[t]*Time);
            }

            // Output the results to List Boxes
            for(int k=0;k<=NK-1;k++)
                for(int t=0;t<NT;t++)
                {
                    Maturities.Items.Add(T[t]);
                    Strikes.Items.Add(K[k]);
                    MktPrices.Items.Add(MktPrice[k,t]);
                    ModelPrices.Items.Add(ModelPrice[k,t]);
                }

            // Output the Market Implied Vol to List View "MktIVOutputList"
            // Output the Model Implied Vol to List View "ModelIVOutputList"
            ListViewItem item;
            for(int k=0;k<=NK-1;k++)
            {
                item = new ListViewItem();
                item.Text = Convert.ToString(Math.Round(MktIV[k,0],4));
                item.SubItems.Add(Convert.ToString(Math.Round(MktIV[k,1],4)));
                item.SubItems.Add(Convert.ToString(Math.Round(MktIV[k,2],4)));
                item.SubItems.Add(Convert.ToString(Math.Round(MktIV[k,3],4)));
                ModelIVOutputList.Items.Add(item);

                item = new ListViewItem();
                item.Text = Convert.ToString(Math.Round(ModelIV[k,0],4));
                item.SubItems.Add(Convert.ToString(Math.Round(ModelIV[k,1],4)));
                item.SubItems.Add(Convert.ToString(Math.Round(ModelIV[k,2],4)));
                item.SubItems.Add(Convert.ToString(Math.Round(ModelIV[k,3],4)));
                MktIVOutputList.Items.Add(item);

                item = new ListViewItem();
                item.Text = Convert.ToString(K[k]);
                StrikeOutputList.Items.Add(item);

            }
        }

        public void WriteToExcel_MouseClick(object sender,MouseEventArgs e)
        {
            Excel.Application xlfile = new Excel.Application();
            Excel.Workbook wb;
            Excel.Worksheet ws;

            xlfile.Visible = false;
            wb = xlfile.Workbooks.Open("C:\\Users\\Fabrice\\Desktop\\Heston Book\\Chapter 2 Model Issues\\C# Code\\Implied_Volatility\\ModelandMarketPrices.xls");
            ws = (Excel.Worksheet)wb.Worksheets[1];

            xlfile.Cells[1,1] = "Market Prices";
            xlfile.Cells[1,7] = "Model Prices";
            xlfile.Cells[1,12] = "Market IV";
            xlfile.Cells[1,17] = "Model IV";
            xlfile.Cells[1,2] = "Maturity";
            xlfile.Cells[2,1] = "Strike";
            for (int t=0; t<=NT-1; t++)
                for(int k=0;k<=NK-1;k++)
                {
                    xlfile.Cells[2,t+2] = T[t];
                    xlfile.Cells[k+3,1] = K[k];
                    xlfile.Cells[k+3,t+2] = MktPrice[k,t];
                    xlfile.Cells[2,t+7] = T[t];
                    xlfile.Cells[k+3,t+7] = ModelPrice[k,t];
                    xlfile.Cells[2,t+12] = T[t];
                    xlfile.Cells[k+3,t+12] = MktIV[k,t];
                    xlfile.Cells[2,t+17] = T[t];
                    xlfile.Cells[k+3,t+17] = ModelIV[k,t];
                }
            wb.Save();
            xlfile.Quit();
        }

        private void ModelIVOutputList_SelectedIndexChanged(object sender,EventArgs e)
        {

        }

    }
}
