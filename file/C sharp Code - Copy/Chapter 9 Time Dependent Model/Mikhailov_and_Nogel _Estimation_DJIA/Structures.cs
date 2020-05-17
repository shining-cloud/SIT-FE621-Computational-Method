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
    public double[] MktIV;     // Implied volatility
    public double[] MktPrice;  // Prices
    public string[] PutCall;   // "P"ut or "C"all 
}
