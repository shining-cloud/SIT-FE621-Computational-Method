// Heston parameters
public struct HParam
{
    public double kappa;        // Mean reversion speed
    public double theta;        // Mean reversion level
    public double sigma;        // Volatility of variance
    public double v0;           // Initial variance
    public double rho;          // Correlation
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

// Settings for the option price calculation
public struct OPSet
{
    public double S;            // Spot price
    public double r;            // Risk free rate
    public double q;            // Dividend
    public int trap;            // 1="Little Trap" characteristic function; 2=Original Heston c.f.
}

// Settings for the objective function
public struct OFSet
{
    public double[] x;
    public double r;
    public double q;
    public double dt;
    public int method;
}
// Settings for the Nelder Mead algorithm
public struct NMSet
{
    public OFSet ofsettings;
    public int MaxIters;
    public double Tolerance;
    public int N; 
}


