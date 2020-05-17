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
    public double[] K;         // Strikes
    public string PutCall;     // "P"ut or "C"all
}

// Settings for the option price calculation
public struct OPSet
{
    public double S;            // Spot price
    public double r;            // Risk free rate
    public double q;            // Dividend
    public int trap;            // 1="Little Trap" characteristic function; 2=Original Heston c.f.
}

// Settings for the objective function with Elices prices
public struct OFSetE
{
    public double[][] paramfixed;       // Fixed parameters from past maturities
    public double v0;                   // Initial variance parameter
    public OPSet opsettings;            // Option settings
    public MktData data;                // Market data
    public double[] T;                  // Maturities
    public double[] X;                  // Gauss Laguerre absicssas
    public double[] W;                  // Gauss Laguerre weights
    public int LossFunction;            // Choice of loss function
    public double[] lb;                 // Lower bounds on parameters
    public double[] ub;                 // Upper bounds on parameters
    public int t;                       // Maturity number
}
// Settings for the objective function with Heston prices
public struct OFSetH
{
    public OPSet opsettings;            // Option settings
    public double[,] MktPrice;          // Market prices
    public double[,] MktIV;             // Market implied vol
    public double[] K;                  // Strikes
    public string PutCall;              // PutCall indicator
    public double[] T;                  // Maturities
    public double[] X;                  // Gauss Laguerre absicssas
    public double[] W;                  // Gauss Laguerre weights
    public int LossFunction;            // Choice of loss function
    public double[] lb;                 // Lower bounds on parameters
    public double[] ub;                 // Upper bounds on parameters
}

// Settings for the Nelder Mead algorithm
public struct NMSet
{
    public OFSetE ofsettingsE;
    public OFSetH ofsettingsH;
    public string choice;
    public int MaxIters;
    public double Tolerance;
    public int N;
}
