 function y = ADIPrice(scheme,thet,params,S0,V0,K,r,q,S,V,T,GridType)

% Alternating Direction Implicit (ADI) scheme for the Heston model.
% INPUTS
%   scheme = 'DO' Douglas, 'CS' Craig-Sneyd, 
%            'MSC' Modified CS, 'HV' Hundsdorfer-Verwer
%   thet  = weighing parameter 0 = Explicit Scheme
%           1 = Implicit Scheme, 0.5 = Crank-Nicolson
%  params = Heston parameters
%  S0 = Spot price on which to price
%  V0 = Volatility on which to price
%  K = Strike price
%  r = Risk free rate
%  q = Dividend yield
%  S = Stock price grid - uniform
%  V = Volatility grid - uniform
%  T = Maturity grid - uniform
%  GridType = 'Uniform' or 'NonUniform'

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);

% Length of grids;
NS = length(S);
NV = length(V);
NT = length(T);
dt = T(2)-T(1);

% Build the matrices for the derivatives and other matrices
% Default setting is Uniform grid
if strcmp(GridType,'NonUniform')
    [derS derSS derV1 derV2 derVV derSV R] = BuildDerivativesNonUniform(S,V,T);
else
    [derS derSS derV1 derV2 derVV derSV R] = BuildDerivatives(S,V,T);
end
    
% Decompose the derivatives matrices
A0 = rho.*sigma.*derSV;
A1 = (r-q).*derS + (1/2).*derSS - r.*R./2;
A2 = kappa.*theta.*derV1 - kappa.*derV2 + (1/2).*sigma^2*derVV - r.*R./2;

% Initialize the u vector, create identity matrix
% u plays the role of U(t-1)
u = zeros(NS*NV,1);
I = eye(NS*NV);

% U(0) vector - value of U(T) at maturity
Si = repmat(S',NV,1);
U = max(0, Si - K);

% Loop through the time increment
for t=2:NT
    u  = U;
    % Vectors common to all schemes
    Y0 = (I + dt.*(A0+A1+A2))*u;
    Y1 = (I - thet.*dt.*A1) \ (Y0 - thet.*dt.*A1*u);
    Y2 = (I - thet.*dt.*A2) \ (Y1 - thet.*dt.*A2*u);
    if strcmp(scheme,'DO')
        % Douglas ADI scheme
        U  = Y2;
    elseif strcmp(scheme,'CS')
        % Craig-Sneyd ADI scheme
        Y0_ = Y0 + (1/2).*dt.*(A0*Y2 - A0*u);
        Y1_ = (I - thet.*dt.*A1) \ (Y0_ - thet.*dt.*A1*u);
        Y2_ = (I - thet.*dt.*A2) \ (Y1_ - thet.*dt.*A2*u);
        U = Y2_;
    elseif strcmp(scheme,'MCS')
        % Modified Craig-Sneyd ADI scheme
        Y0h = Y0 + thet.*dt.*(A0*Y2 - A0*u);
        Y0_ = Y0h + (1/2-thet).*dt.*((A0+A1+A2)*Y2 - (A0+A1+A2)*u);
        Y1_ = (I - thet.*dt.*A1) \ (Y0_ - thet.*dt.*A1*u);
        Y2_ = (I - thet.*dt.*A2) \ (Y1_ - thet.*dt.*A2*u);
        U = Y2_;
    elseif strcmp(scheme,'HV')
        % Hundsdorfer-Verwer ADI scheme
        Y0_ = Y0 + (1/2).*dt.*((A0+A1+A2)*Y2 - (A0+A1+A2)*u);
        Y1_ = (I - thet.*dt.*A1) \ (Y0_ - thet.*dt.*A1*Y2);
        Y2_ = (I - thet.*dt.*A2) \ (Y1_ - thet.*dt.*A2*Y2);
        U = Y2_;
    end
    disp(['Completed Time Iteration Step ' num2str(t) ' of ' num2str(NT) ])
end

% Restack the U vector to output a matrix
U = reshape(U,NS,NV);

% Interpolate to get the price at S0 and v0
y = interp2(V,S,U,V0,S0);
