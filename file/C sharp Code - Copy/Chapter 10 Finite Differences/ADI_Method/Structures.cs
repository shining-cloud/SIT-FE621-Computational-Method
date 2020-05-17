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

// Sub matrices of the "L" matrix
public struct LMatrices
{
    public double[,] derS;
    public double[,] derSS;
    public double[,] derV1;
    public double[,] derV2;
    public double[,] derVV;
    public double[,] derSV;
    public double[,] R;
}

// Lower and Upper diagonal matrices
public struct LUstruct
{
    public double[,] LM;
    public double[,] UM;
}
