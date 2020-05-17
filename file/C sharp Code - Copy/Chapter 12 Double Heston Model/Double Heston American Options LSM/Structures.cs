// Simulation results
public struct DHSim
{
    public double[,] S;         // Stock price
    public double EuroPrice;
}
// Double Heston parameters
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
// Single Heston parameters
public struct HParam
{
    public double kappa;        // Mean reversion speed
    public double theta;        // Mean reversion level
    public double sigma;        // Volatility of variance
    public double v0;           // Initial variance
    public double rho;          // Correlation
}
// Option settings
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
// LU decomposition
public struct LUstruct
{
    public double[,] LM;
    public double[,] UM;
}
