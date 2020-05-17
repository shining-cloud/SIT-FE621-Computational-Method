// Bivariate tree
struct BivariateStruct
{
    public float[,] Euro;
    public float[,] Amer;
    public float[,] Yt;
    public float[,] V;
    public float[,] X;
    public float[,] Prob;
    public float EuroPrice;
    public float AmerPrice;
}

// Heston parameters
public struct HParam
{
    public float kappa;        // Mean reversion speed
    public float theta;        // Mean reversion level
    public float sigma;        // Volatility of variance
    public float v0;           // Initial variance
    public float rho;          // Correlation
}
