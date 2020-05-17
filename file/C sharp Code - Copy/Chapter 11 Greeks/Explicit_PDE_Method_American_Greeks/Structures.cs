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
// Return U(S,v,T) and U(S,v,T-dt)
public struct Uu
{
    public double[,] bigU;
    public double[,] smallU;
}
