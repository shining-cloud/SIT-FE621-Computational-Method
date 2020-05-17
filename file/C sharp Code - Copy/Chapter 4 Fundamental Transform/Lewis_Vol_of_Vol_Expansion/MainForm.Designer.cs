namespace Lewis_Vol_of_Vol_Expansion
{
    partial class MainForm
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
            this.GenerateLewisTable = new System.Windows.Forms.Button();
            this.Strike = new System.Windows.Forms.ListBox();
            this.ExactPrice = new System.Windows.Forms.ListBox();
            this.IVExact = new System.Windows.Forms.ListBox();
            this.SeriesICall = new System.Windows.Forms.ListBox();
            this.IV1 = new System.Windows.Forms.ListBox();
            this.SeriesIICall = new System.Windows.Forms.ListBox();
            this.IV2 = new System.Windows.Forms.ListBox();
            this.BlackScholesPrice = new System.Windows.Forms.ListBox();
            this.IVBS = new System.Windows.Forms.ListBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label7 = new System.Windows.Forms.Label();
            this.label8 = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // GenerateLewisTable
            // 
            this.GenerateLewisTable.Location = new System.Drawing.Point(436,26);
            this.GenerateLewisTable.Name = "GenerateLewisTable";
            this.GenerateLewisTable.Size = new System.Drawing.Size(286,40);
            this.GenerateLewisTable.TabIndex = 0;
            this.GenerateLewisTable.Text = "Generate Table 3.3.1 of Lewis (2001)";
            this.GenerateLewisTable.UseVisualStyleBackColor = true;
            this.GenerateLewisTable.MouseClick += new System.Windows.Forms.MouseEventHandler(this.GenerateLewisTable_MouseClick);
            // 
            // Strike
            // 
            this.Strike.FormattingEnabled = true;
            this.Strike.Location = new System.Drawing.Point(30,122);
            this.Strike.Name = "Strike";
            this.Strike.Size = new System.Drawing.Size(106,186);
            this.Strike.TabIndex = 1;
            // 
            // ExactPrice
            // 
            this.ExactPrice.FormattingEnabled = true;
            this.ExactPrice.Location = new System.Drawing.Point(149,122);
            this.ExactPrice.Name = "ExactPrice";
            this.ExactPrice.Size = new System.Drawing.Size(106,186);
            this.ExactPrice.TabIndex = 2;
            // 
            // IVExact
            // 
            this.IVExact.FormattingEnabled = true;
            this.IVExact.Location = new System.Drawing.Point(270,122);
            this.IVExact.Name = "IVExact";
            this.IVExact.Size = new System.Drawing.Size(106,186);
            this.IVExact.TabIndex = 3;
            // 
            // SeriesICall
            // 
            this.SeriesICall.FormattingEnabled = true;
            this.SeriesICall.Location = new System.Drawing.Point(391,122);
            this.SeriesICall.Name = "SeriesICall";
            this.SeriesICall.Size = new System.Drawing.Size(106,186);
            this.SeriesICall.TabIndex = 4;
            // 
            // IV1
            // 
            this.IV1.FormattingEnabled = true;
            this.IV1.Location = new System.Drawing.Point(515,122);
            this.IV1.Name = "IV1";
            this.IV1.Size = new System.Drawing.Size(106,186);
            this.IV1.TabIndex = 5;
            // 
            // SeriesIICall
            // 
            this.SeriesIICall.FormattingEnabled = true;
            this.SeriesIICall.Location = new System.Drawing.Point(639,122);
            this.SeriesIICall.Name = "SeriesIICall";
            this.SeriesIICall.Size = new System.Drawing.Size(106,186);
            this.SeriesIICall.TabIndex = 6;
            // 
            // IV2
            // 
            this.IV2.FormattingEnabled = true;
            this.IV2.Location = new System.Drawing.Point(762,122);
            this.IV2.Name = "IV2";
            this.IV2.Size = new System.Drawing.Size(106,186);
            this.IV2.TabIndex = 7;
            // 
            // BlackScholesPrice
            // 
            this.BlackScholesPrice.FormattingEnabled = true;
            this.BlackScholesPrice.Location = new System.Drawing.Point(886,122);
            this.BlackScholesPrice.Name = "BlackScholesPrice";
            this.BlackScholesPrice.Size = new System.Drawing.Size(106,186);
            this.BlackScholesPrice.TabIndex = 8;
            // 
            // IVBS
            // 
            this.IVBS.FormattingEnabled = true;
            this.IVBS.Location = new System.Drawing.Point(1013,122);
            this.IVBS.Name = "IVBS";
            this.IVBS.Size = new System.Drawing.Size(106,186);
            this.IVBS.TabIndex = 9;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(405,95);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(62,13);
            this.label1.TabIndex = 10;
            this.label1.Text = "Series I Call";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(525,95);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(96,13);
            this.label2.TabIndex = 11;
            this.label2.Text = "Series I Implied Vol";
            this.label2.Click += new System.EventHandler(this.label2_Click);
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(657,95);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(65,13);
            this.label3.TabIndex = 12;
            this.label3.Text = "Series II Call";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(772,95);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(96,13);
            this.label4.TabIndex = 13;
            this.label4.Text = "Series I Implied Vol";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(890,95);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(102,13);
            this.label5.TabIndex = 14;
            this.label5.Text = "Black Scholes Price";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(52,95);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(34,13);
            this.label6.TabIndex = 15;
            this.label6.Text = "Strike";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(177,95);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(61,13);
            this.label7.TabIndex = 16;
            this.label7.Text = "Exact Price";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(288,95);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(88,13);
            this.label8.TabIndex = 17;
            this.label8.Text = "Exact Implied Vol";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(1010,95);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(129,13);
            this.label9.TabIndex = 18;
            this.label9.Text = "Black Scholes Implied Vol";
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F,13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1153,392);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.IVBS);
            this.Controls.Add(this.BlackScholesPrice);
            this.Controls.Add(this.IV2);
            this.Controls.Add(this.SeriesIICall);
            this.Controls.Add(this.IV1);
            this.Controls.Add(this.SeriesICall);
            this.Controls.Add(this.IVExact);
            this.Controls.Add(this.ExactPrice);
            this.Controls.Add(this.Strike);
            this.Controls.Add(this.GenerateLewisTable);
            this.Name = "MainForm";
            this.Text = "Form1";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button GenerateLewisTable;
        private System.Windows.Forms.ListBox Strike;
        private System.Windows.Forms.ListBox ExactPrice;
        private System.Windows.Forms.ListBox IVExact;
        private System.Windows.Forms.ListBox SeriesICall;
        private System.Windows.Forms.ListBox IV1;
        private System.Windows.Forms.ListBox SeriesIICall;
        private System.Windows.Forms.ListBox IV2;
        private System.Windows.Forms.ListBox BlackScholesPrice;
        private System.Windows.Forms.ListBox IVBS;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label label9;
    }
}

