namespace Implied_Volatility
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if(disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.GeneratePrices = new System.Windows.Forms.Button();
            this.MktPrices = new System.Windows.Forms.ListBox();
            this.ModelPrices = new System.Windows.Forms.ListBox();
            this.Strikes = new System.Windows.Forms.ListBox();
            this.Maturities = new System.Windows.Forms.ListBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.WriteToExcel = new System.Windows.Forms.Button();
            this.ModelIVOutputList = new System.Windows.Forms.ListView();
            this.Mat1 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.Mat2 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.Mat3 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.Mat4 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.label5 = new System.Windows.Forms.Label();
            this.MktIVOutputList = new System.Windows.Forms.ListView();
            this.Maturity1 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.Maturity2 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.Maturity3 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.Maturity4 = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.label6 = new System.Windows.Forms.Label();
            this.StrikeOutputList = new System.Windows.Forms.ListView();
            this.Strike = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.SuspendLayout();
            // 
            // GeneratePrices
            // 
            this.GeneratePrices.Location = new System.Drawing.Point(28,27);
            this.GeneratePrices.Name = "GeneratePrices";
            this.GeneratePrices.Size = new System.Drawing.Size(178,27);
            this.GeneratePrices.TabIndex = 0;
            this.GeneratePrices.Text = "Generate Prices";
            this.GeneratePrices.UseVisualStyleBackColor = true;
            this.GeneratePrices.MouseClick += new System.Windows.Forms.MouseEventHandler(this.GeneratePrices_MouseClick_1);
            // 
            // MktPrices
            // 
            this.MktPrices.FormattingEnabled = true;
            this.MktPrices.Location = new System.Drawing.Point(309,107);
            this.MktPrices.Name = "MktPrices";
            this.MktPrices.Size = new System.Drawing.Size(120,407);
            this.MktPrices.TabIndex = 1;
            // 
            // ModelPrices
            // 
            this.ModelPrices.FormattingEnabled = true;
            this.ModelPrices.Location = new System.Drawing.Point(457,107);
            this.ModelPrices.Name = "ModelPrices";
            this.ModelPrices.Size = new System.Drawing.Size(121,407);
            this.ModelPrices.TabIndex = 3;
            // 
            // Strikes
            // 
            this.Strikes.FormattingEnabled = true;
            this.Strikes.Location = new System.Drawing.Point(162,107);
            this.Strikes.Name = "Strikes";
            this.Strikes.Size = new System.Drawing.Size(120,407);
            this.Strikes.TabIndex = 4;
            // 
            // Maturities
            // 
            this.Maturities.FormattingEnabled = true;
            this.Maturities.Location = new System.Drawing.Point(28,107);
            this.Maturities.Name = "Maturities";
            this.Maturities.Size = new System.Drawing.Size(120,407);
            this.Maturities.TabIndex = 5;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(318,83);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(72,13);
            this.label1.TabIndex = 6;
            this.label1.Text = "Market Prices";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(486,83);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(68,13);
            this.label2.TabIndex = 7;
            this.label2.Text = "Model Prices";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(190,83);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(39,13);
            this.label3.TabIndex = 8;
            this.label3.Text = "Strikes";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(68,83);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(52,13);
            this.label4.TabIndex = 9;
            this.label4.Text = "Maturities";
            // 
            // WriteToExcel
            // 
            this.WriteToExcel.Location = new System.Drawing.Point(457,27);
            this.WriteToExcel.Name = "WriteToExcel";
            this.WriteToExcel.Size = new System.Drawing.Size(121,23);
            this.WriteToExcel.TabIndex = 10;
            this.WriteToExcel.Text = "Write to Excel";
            this.WriteToExcel.UseVisualStyleBackColor = true;
            this.WriteToExcel.MouseClick += new System.Windows.Forms.MouseEventHandler(this.WriteToExcel_MouseClick);
            // 
            // ModelIVOutputList
            // 
            this.ModelIVOutputList.Columns.AddRange(new System.Windows.Forms.ColumnHeader[] {
            this.Mat1,
            this.Mat2,
            this.Mat3,
            this.Mat4});
            this.ModelIVOutputList.Location = new System.Drawing.Point(628,107);
            this.ModelIVOutputList.Name = "ModelIVOutputList";
            this.ModelIVOutputList.Size = new System.Drawing.Size(248,496);
            this.ModelIVOutputList.TabIndex = 11;
            this.ModelIVOutputList.UseCompatibleStateImageBehavior = false;
            this.ModelIVOutputList.View = System.Windows.Forms.View.Details;
            this.ModelIVOutputList.SelectedIndexChanged += new System.EventHandler(this.ModelIVOutputList_SelectedIndexChanged);
            // 
            // Mat1
            // 
            this.Mat1.Text = "Maturity 1";
            // 
            // Mat2
            // 
            this.Mat2.Text = "Maturity 2";
            // 
            // Mat3
            // 
            this.Mat3.Text = "Maturity 3";
            // 
            // Mat4
            // 
            this.Mat4.Text = "Maturity 4";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(625,82);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(90,13);
            this.label5.TabIndex = 12;
            this.label5.Text = "Model Implied Vol";
            // 
            // MktIVOutputList
            // 
            this.MktIVOutputList.Columns.AddRange(new System.Windows.Forms.ColumnHeader[] {
            this.Maturity1,
            this.Maturity2,
            this.Maturity3,
            this.Maturity4});
            this.MktIVOutputList.Location = new System.Drawing.Point(972,107);
            this.MktIVOutputList.Name = "MktIVOutputList";
            this.MktIVOutputList.Size = new System.Drawing.Size(245,496);
            this.MktIVOutputList.TabIndex = 13;
            this.MktIVOutputList.UseCompatibleStateImageBehavior = false;
            this.MktIVOutputList.View = System.Windows.Forms.View.Details;
            // 
            // Maturity1
            // 
            this.Maturity1.Text = "Maturity 1";
            // 
            // Maturity2
            // 
            this.Maturity2.Text = "Maturity 2";
            // 
            // Maturity3
            // 
            this.Maturity3.Text = "Maturity 3";
            // 
            // Maturity4
            // 
            this.Maturity4.Text = "Maturity 4";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(969,82);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(94,13);
            this.label6.TabIndex = 14;
            this.label6.Text = "Market Implied Vol";
            // 
            // StrikeOutputList
            // 
            this.StrikeOutputList.Columns.AddRange(new System.Windows.Forms.ColumnHeader[] {
            this.Strike});
            this.StrikeOutputList.Location = new System.Drawing.Point(894,107);
            this.StrikeOutputList.Name = "StrikeOutputList";
            this.StrikeOutputList.Size = new System.Drawing.Size(65,496);
            this.StrikeOutputList.TabIndex = 15;
            this.StrikeOutputList.UseCompatibleStateImageBehavior = false;
            this.StrikeOutputList.View = System.Windows.Forms.View.Details;
            // 
            // Strike
            // 
            this.Strike.Text = "Strike";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F,13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1256,632);
            this.Controls.Add(this.StrikeOutputList);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.MktIVOutputList);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.ModelIVOutputList);
            this.Controls.Add(this.WriteToExcel);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.Maturities);
            this.Controls.Add(this.Strikes);
            this.Controls.Add(this.ModelPrices);
            this.Controls.Add(this.MktPrices);
            this.Controls.Add(this.GeneratePrices);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button GeneratePrices;
        private System.Windows.Forms.ListBox MktPrices;
        private System.Windows.Forms.ListBox ModelPrices;
        private System.Windows.Forms.ListBox Strikes;
        private System.Windows.Forms.ListBox Maturities;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Button WriteToExcel;
        private System.Windows.Forms.ListView ModelIVOutputList;
        private System.Windows.Forms.ColumnHeader Mat1;
        private System.Windows.Forms.ColumnHeader Mat2;
        private System.Windows.Forms.ColumnHeader Mat3;
        private System.Windows.Forms.ColumnHeader Mat4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.ListView MktIVOutputList;
        private System.Windows.Forms.ColumnHeader Maturity1;
        private System.Windows.Forms.ColumnHeader Maturity2;
        private System.Windows.Forms.ColumnHeader Maturity3;
        private System.Windows.Forms.ColumnHeader Maturity4;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.ListView StrikeOutputList;
        private System.Windows.Forms.ColumnHeader Strike;
    }
}

