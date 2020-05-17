// Heston parameters
public struct HParam
{
    public double kappa;        // Mean reversion speed
    public double theta;        // Mean reversion level
    public double sigma;        // Volatility of variance
    public double v0;           // Initial variance
    public double rho;          // Correlation
    public double lambda;       // Risk parameter
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

// Quantities returned by the multi-domain Gauss-Laguerre integration
public struct OutputMD
{
    public double Price;       // The option price
    public double lower;       // The lower limit of the integration domain
    public double upper;       // The upper limit of the integration domain
    public int Npoints;        // Number of integration points
}

