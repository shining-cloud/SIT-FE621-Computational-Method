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

// Settings for the option price 
public struct OpSet
{
    public double S;            // Spot price
    public double K;            // Strike price
    public double T;            // Maturity
    public double r;            // Risk free rate
    public double q;            // Dividend
    public string PutCall;      // "P"ut or "C"all
    public int trap;            // 1="Little Trap" characteristic function; 2=Original Heston c.f.
}

// Settings for the Medvedev-Scaillet expansion
public struct MSSet
{
    public int method;          // Newton Cotes integration method
    public double A;            // Newton Cotes lower limit
    public double B;            // Newton Cotes upper limit
    public int N;               // Newton Cotes number of points
    public double dt;           // Bisection method increment for derivative
    public double tol;          // Bisection method tolerance
    public int MaxIter;         // Bisection maximum number of interations
    public int NumTerms;        // M.S. number of terms in expansion
    public double yinf;         // M.S. y infinite for European put
    public double a;            // M.S. lower bound for bisection method to find y (-10)
    public double b;            // M.S. upper bound for bisection method to find y (10)
}

// MS Black Scholes
public struct MSsetBS
{
    public double theta;
    public double K;
    public double sigma;
    public double r;
    public double q;
    public double T;
}
// Settings for the Market data;
public struct MktData
{
    public double[] MktIV;     // Implied volatility
    public double[] K;         // Strikes
    public string PutCall;     // "P"ut or "C"all
}

// Settings for the objective function
public struct OFSet
{
//    public HParam param;
    public OpSet opsettings;
    public MSSet mssettings;
    public MktData data;
    public double[] X;
    public double[] W;
    public double[] lb;
    public double[] ub;
}
// Settings for the Nelder Mead algorithm
public struct NMSet
{
    public OFSet ofsettings;
    public int MaxIters;
    public double Tolerance;
    public int N;
}
