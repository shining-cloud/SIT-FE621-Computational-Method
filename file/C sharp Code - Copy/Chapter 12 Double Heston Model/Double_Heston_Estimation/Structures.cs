// Heston parameters
public struct DHParam
{
    public double kappa1;        // Mean reversion speed
    public double theta1;        // Mean reversion level
    public double sigma1;        // Volatility of variance
    public double v01;           // Initial variance
    public double rho1;          // Correlation
    public double kappa2;        // Mean reversion speed
    public double theta2;        // Mean reversion level
    public double sigma2;        // Volatility of variance
    public double v02;           // Initial variance
    public double rho2;          // Correlation
}

// Settings for the option price calculation
public struct OpSet
{
    public double S;            // Spot price
    public double K;            // Strke price
    public double T;            // Maturity
    public double r;            // Risk free rate
    public double q;            // Dividend
    public string PutCall;      // "P"ut or "C"all
    public int trap;            // 1="Little Trap" characteristic function; 2=Original Heston c.f.
}

// Settings for the Market data;
public struct MktData
{
    public double[,] MktIV;     // Implied volatility
    public double[,] MktPrice;  // Prices
    public double[] K;          // Strikes
    public double[] T;          // Maturities
    public string[,] PutCall;   // "P"ut or "C"all
}
