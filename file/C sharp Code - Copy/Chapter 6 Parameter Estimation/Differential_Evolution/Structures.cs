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
public struct OPSet
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

// Settings for the objective function
public struct OFSet
{
    public OPSet opset;
    public MktData data;
    public double[] X;
    public double[] W;
    public int LossFunction;
    public double[] lb;
    public double[] ub;
    public string CF;
}

// Settings for the Nelder Mead algorithm
public struct NMSet
{
    public OFSet ofset;
    public int MaxIters;
    public double Tolerance;
    public int N;
}

// Settings for differential evolution
public struct DEParam
{
    public int NG;              // Number of generations (iterations)
    public int NP;              // Number of population members
    public double CR;           // Crossover ratio (=0.5)
    public double F;            // Threshold (=0.8)
    public double[] ub;         // Vector of upper bounds for the Heston parameters
    public double[] lb;         // Vector of lower bounds for the parameters
}


