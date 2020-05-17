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
// B0 and B1 vectors for Chiarella
public struct BVec
{
    public double[] B0;
    public double[] B1;
}
// Settings for the option
public struct OPSet
{
    public double K;
    public double S;
    public double r;
    public double q;
}
