// Heston parameters
public struct HParam
{
    public double kappa;        // Mean reversion speed
    public double theta;        // Mean reversion level
    public double sigma;        // Volatility of variance
    public double v0;           // Initial variance
    public double rho;          // Correlation
    public double lambda;       // Risk
}
// Settings for the option price calculation
public struct OpSet
{
    public double S;            // Spot price
    public double r;            // Risk free rate
    public double q;            // Dividend
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
