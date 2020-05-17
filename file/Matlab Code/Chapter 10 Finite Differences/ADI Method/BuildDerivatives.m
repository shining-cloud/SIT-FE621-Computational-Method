function [derS derSS derV1 derV2 derVV derSV R] = BuildDerivatives(S,V,T)

% Build the first- and second-order derivatives for the "L" operator matrix
% for the Weighted method
% Requires a uniform grid for S and V
% INPUTS
%    S = vector for uniform stock price grid
%    V = vector for uniform volatility grid
%    T = vector for uniform maturity grid
%    thet = parameter for the weighted scheme
% OUTPUTS
%    Matrices of dimnension N x N (N=NS+NV)
%    derS  = Matrix for first-order derivative dU/dS 
%    derSS = Matrix for second-order derivative dU2/dS2 
%    derV1 = Matrix for first-order derivative dU/dV, kappa*theta portion
%    derV2 = Matrix for first-order derivative dU/dV, -kappa*V portion
%    derVV = Matrix for first-order derivative dU2/dV2
%    derSV = Matrix for first-order derivative dU2/dSdV
%    R     = Matrix for r*S(v,t) portion of the PDE

% Length of stock price, volatility, and maturity
NS = length(S);
NV = length(V);
NT = length(T);
Smin = S(1);  Smax = S(NS);
Vmin = V(1);  Vmax = V(NV);
Tmin = T(1);  Tmax = T(NT);

% Increment for Stock Price, Volatility, and Maturity
ds = (Smax-Smin)/(NS-1);
dv = (Vmax-Vmin)/(NV-1);
dt = (Tmax-Tmin)/(NT-1);

% Preliminary quantities
% Size of the U(t) vector and L matrix
N = NS*NV;

% The vectors for S and V, stacked
Si = repmat(S',NV,1);
Vi = reshape(kron(V,ones(NS,1)),N,1);

% Identification of the boundary points
VminB = zeros(N,1);
VmaxB = zeros(N,1);
SminB = zeros(N,1);
SmaxB = zeros(N,1);

% Vmin and Vmax
VminB(1:NS-1) = 1; 
VmaxB(N-NS+2:N) = 1;

% Smin and Smax
for i=NS+1:N
	if mod(i,NS)==0 & (i~=N)
		SmaxB(i) = 1;        % Indices for Smax
	end
	if mod(i,NS)==1  %& (i~=1) 
		SminB(i) = 1;        % Indices for Smin
	end
end

% Identification of the non-boundary points
NB = zeros(N,1);
for b=2:NV-1
	for k=b*NS-(NS-2):b*NS-1;
		NB(k) = 1;
	end
end
NB(NS) = 1;

% Indices for S
Cs = zeros(N,1);
Fs = zeros(N,1);
Bs = zeros(N,1);

for b=2:NV-1
	for k=b*NS-(NS-3):b*NS-2;
		Cs(k) = 1;           % Central differences for S
	end
end
Fs((2:NV-1)*NS-(NS-2)) = 1;  % Backward differences for S
Bs((2:NV-1)*NS-1) = 1;       % Forward differences for S

% Indices for V
Cv = zeros(N,1);
Fv = zeros(N,1);
Bv = zeros(N,1);

for b=3:NV-2
	for k=b*NS-(NS-2):b*NS-1;
        Cv(k) = 1;                          % Central differences for V
	end
end
Fv(2*NS-(NS-2):2*NS-1) = 1;                 % Forward differences for V
Bv((NV-1)*NS-(NS-2):(NV-1)*NS-1) = 1;       % Backward differences for V
  
% Indices for SV
Csv = zeros(N,1);
for b=2:NV-1
	for k=b*NS-(NS-2):b*NS-1
		Csv(k) = 1;                         % Central difference for SV
	end
end

%[Si Vi SminB SmaxB VminB VmaxB NB Fs Bs Cs Fv Bv Cv Csv]

% Create the matrices for the derivatives
derS  = zeros(N,N);
derSS = zeros(N,N);
derV1 = zeros(N,N);
derV2 = zeros(N,N);
derVV = zeros(N,N);
derSV = zeros(N,N);
R     = eye(N);

% NON BOUNDARY POINTS ----------------------------------
I = find(Cs==1);
for k=1:length(I)
	% Create the matrix for S-derivatives                   % Central differences
	derS(I(k),I(k)-1)  =  -1/2/ds * Si(I(k));               % U(s-1,v)
	derS(I(k),I(k))    =   0;                               % U(s,v)
	derS(I(k),I(k)+1)  =   1/2/ds * Si(I(k));               % U(s+1,v)
	derSS(I(k),I(k)-1) =   1/ds^2 * Vi(I(k))*Si(I(k))^2;    % U(s-1,v)
	derSS(I(k),I(k))   =  -2/ds^2 * Vi(I(k))*Si(I(k))^2;    % U(s,v)
	derSS(I(k),I(k)+1) =   1/ds^2 * Vi(I(k))*Si(I(k))^2;    % U(s+1,v)
end
I = find(Fs==1);
for k=1:length(I)                                           % Forward differences
	derS(I(k),I(k))    = -3/2/ds  * Si(I(k));               % U(s,v)
	derS(I(k),I(k)+1)  =    2/ds  * Si(I(k));               % U(s+1,v)
	derS(I(k),I(k)+2)  =   -1/ds  * Si(I(k));               % U(s+2,v)
	derSS(I(k),I(k))   =    1/ds^2 * Vi(I(k))*Si(I(k))^2;   % U(s,v)
	derSS(I(k),I(k)+1) =   -2/ds^2 * Vi(I(k))*Si(I(k))^2;   % U(s+1,v)
	derSS(I(k),I(k)+2) =    1/ds^2 * Vi(I(k))*Si(I(k))^2;   % U(s+2,v)
end
I = find(Bs==1);
for k=1:length(I)                                           % Backward differences
	derS(I(k),I(k)-2)  =  -1/2/ds * Si(I(k));               % U(s-2,v)
	derS(I(k),I(k)-1)  =    -2/ds * Si(I(k));               % U(s-1,v)
	derS(I(k),I(k))    =   3/2/ds * Si(I(k));               % U(s,v)
	derSS(I(k),I(k)-2) =   1/ds^2 * Vi(I(k))*Si(I(k))^2;    % U(s-2,v)
	derSS(I(k),I(k)-1) =  -2/ds^2 * Vi(I(k))*Si(I(k))^2;    % U(s-1,v)
	derSS(I(k),I(k))   =   1/ds^2 * Vi(I(k))*Si(I(k))^2;    % U(s,v)
end
% Create the matrix for V-derivatives
I = find(Cv==1);
for k=1:length(I)                                     % Central differences
	derV1(I(k),I(k)-NS)  = -1/2/dv;                   % U(s,v-1)
	derV1(I(k),I(k))     =  0;                        % U(s,v)
	derV1(I(k),I(k)+NS)  =  1/2/dv;                   % U(s,v+1)
	derV2(I(k),I(k)-NS)  = -1/2/dv * Vi(I(k));        % U(s,v-1)
	derV2(I(k),I(k))     =  0;                        % U(s,v)
	derV2(I(k),I(k)+NS)  =  1/2/dv * Vi(I(k));        % U(s,v+1)
	derVV(I(k),I(k)-NS)  =  1/dv^2 * Vi(I(k));        % U(s,v-1)
	derVV(I(k),I(k))     = -2/dv^2 * Vi(I(k));        % U(s,v)
	derVV(I(k),I(k)+NS)  =  1/dv^2 * Vi(I(k));        % U(s,v+1)
end
I = find(Fv==1);
for k=1:length(I)                                     % Forward differences
	derV1(I(k),I(k))       = -3/2/dv;                 % U(s,v)
	derV1(I(k),I(k)+NS)    =    2/dv;                 % U(s,v+1)
	derV1(I(k),I(k)+2*NS)  =   -1/dv;                 % U(s,v+2)
	derV2(I(k),I(k))       = -3/2/dv * Vi(I(k));      % U(s,v)
	derV2(I(k),I(k)+NS)    =    2/dv * Vi(I(k));      % U(s,v+1)
	derV2(I(k),I(k)+2*NS)  =   -1/dv * Vi(I(k));      % U(s,v+2)
	derVV(I(k),I(k))       =    1/dv^2 * Vi(I(k));    % U(s,v)
	derVV(I(k),I(k)+NS)    =   -2/dv^2 * Vi(I(k));    % U(s,v+1)
	derVV(I(k),I(k)+2*NS)  =    1/dv^2 * Vi(I(k));    % U(s,v+2)
end
I = find(Bv==1);
for k=1:length(I);                                    % Backward differences
	derV1(I(k),I(k)-2*NS)  = -1/2/dv;                 % U(s,v-2)
	derV1(I(k),I(k)-NS)    =   -2/dv;                 % U(s,v-1)
	derV1(I(k),I(k))       =  3/2/dv;                 % U(s,v)
	derV2(I(k),I(k)-2*NS)  = -1/2/dv * Vi(I(k));      % U(s,v-2)
	derV2(I(k),I(k)-NS)    =   -2/dv * Vi(I(k));      % U(s,v-1)
	derV2(I(k),I(k))       =  3/2/dv * Vi(I(k));      % U(s,v)
	derVV(I(k),I(k)-2*NS)  =  1/dv^2 * Vi(I(k));      % U(s,v-2)
	derVV(I(k),I(k)-NS)    = -2/dv^2 * Vi(I(k));      % U(s,v-1)
	derVV(I(k),I(k))       =  1/dv^2 * Vi(I(k));      % U(s,v)
end
% Create the matrix for SV-derivatives
I = find(Csv==1);
for k=1:length(I)
	derSV(I(k),I(k)+NS+1) =  1/(4*ds*dv) * Vi(I(k))*Si(I(k));  % U(s+1,v+1)
	derSV(I(k),I(k)+NS-1) = -1/(4*ds*dv) * Vi(I(k))*Si(I(k));  % U(s-1,v+1)
	derSV(I(k),I(k)-NS-1) =  1/(4*ds*dv) * Vi(I(k))*Si(I(k));  % U(s-1,v-1)
	derSV(I(k),I(k)-NS+1) = -1/(4*ds*dv) * Vi(I(k))*Si(I(k));  % U(s+1,v-1)
end

% BOUNDARY POINTS ----------------------------------
% Boundary for Smin
I = find(SminB==1);
for k=1:length(I);
	derS(I(k),I(k))  = 0;
	derSS(I(k),I(k)) = 0;
	derV1(I(k),I(k)) = 0;
	derV2(I(k),I(k)) = 0;
	derVV(I(k),I(k)) = 0;
	derSV(I(k),I(k)) = 0;
	R(I(k),I(k)) = 0;
end
% Boundary condition for Smax
I = find(SmaxB==1);
for k=1:length(I);
	derS(I(k),I(k))     = Si(I(k));
	derSS(I(k),I(k))    = 0;
	derSV(I(k),I(k))    = 0;                     % Central difference
	derV1(I(k),I(k)-NS) = -1/2/dv;               % U(s,v-1)
	derV1(I(k),I(k))    =  0;                    % U(s,v)
	derV1(I(k),I(k)+NS) =  1/2/dv;               % U(s,v+1)
	derV2(I(k),I(k)-NS) = -1/2/dv * Vi(I(k));    % U(s,v-1)
	derV2(I(k),I(k))    =  0;                    % U(s,v)
	derV2(I(k),I(k)+NS) =  1/2/dv * Vi(I(k));    % U(s,v+1)
	derVV(I(k),I(k)-NS) =  1/dv^2 * Vi(I(k));    % U(s,v-1)
	derVV(I(k),I(k))    = -2/dv^2 * Vi(I(k));    % U(s,v)
	derVV(I(k),I(k)+NS) =  1/dv^2 * Vi(I(k));    % U(s,v+1)
end
% Boundary condition for Vmax.  Only the LS submatrix is non-zero
I = find(VmaxB==1);
for k=1:length(I);
	derS(I(k),I(k)) = Si(I(k));
end
% Boundary condition for Vmin
% First row
derS(1,1)        = -3/2/ds  * Si(1);     % U(s,v)
derS(1,2)        =    2/ds  * Si(1);     % U(s+1,v)
derS(1,3)        =   -1/ds  * Si(1);     % U(s+2,v)
derV1(1,1)       = -3/2/dv;              % U(s,v)
derV1(1,1+NS)    =    2/dv;              % U(s,v+1)
derV1(1,1+2*NS)  =   -1/dv;              % U(s,v+2)
% Last row
derS(N-1,N-3)       = -1/2/ds * Si(N-1);      % U(s-2,v)
derS(N-1,N-2)       =   -2/ds * Si(N-1);      % U(s-1,v)
derS(N-1,N-1)       =  3/2/ds * Si(N-1);      % U(s,v)
derV1(NS-1,NS-1+2*NS) = -1/2/dv * Vi(NS-1);      % U(s,v+2)
derV1(NS-1,NS-1+NS)   =    2/dv * Vi(NS-1);      % U(s,v+1)
derV1(NS-1,NS-1)      = -3/2/dv * Vi(NS-1);      % U(s,v)

% Other rows
for i=2:NS-2
	derS(i,i-1)     = -1/2/ds * Si(i);      % U(s-1,v)
	derS(i,i)       =  0;                   % U(s,v)
	derS(i,i+1)     =  1/2/ds * Si(i);      % U(s+1,v)
	derV1(i,i)      = -3/2/dv;              % U(s,v)
	derV1(i,i+NS)   =    2/dv;              % U(s,v+1)
	derV1(i,i+2*NS) =   -1/dv;              % U(s,v+2)
end
