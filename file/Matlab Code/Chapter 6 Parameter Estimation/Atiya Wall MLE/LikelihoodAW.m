function [y v] = LikelihoodAW(param,x,r,q,dt,method)

% INPUTS
%   param = vector of parameters
%   x = vector of log-stock prices
%   r = risk free rate
%   q = dividend yield
%   dt = time increment (example: 1/250 is daily)
%   method
%    1 = Likelihood
%    2 = log-likelihood
% OUTPUT
%   the negative likelihood or negative log-likelihood

% Name the Heston parameters
kappa = param(1);
theta = param(2);
sigma = param(3);
v0    = param(4);
rho   = param(5);

% Atiya and Wall parameterization
alpha = kappa*theta;
beta  = kappa;

% Number of log-stock prices
T = length(x);

% Drift term
mu = r - q;

% Equation (17)
betap = 1 - beta*dt;

% Equation (18) - denominator of d(t)
D = 2*pi*sigma*sqrt(1-rho^2)*dt;

% Equation (14)
a = (betap^2 + rho*sigma*betap*dt + sigma^2*dt^2/4) / (2*sigma^2*(1-rho^2)*dt);

% Variance and likelihood at time t = 0
v(1) = v0;
if method==1
	L(1) = exp(-v(1));    % Construct the Likelihood
elseif method==2
	L(1) = -v(1);         % Construct the log-likelihood
end

% Construction the likelihood for time t = 1 through t = T
for t=1:T-1
	% Stock price increment
	dx  = x(t+1) - x(t);
	% Equations (31) and (32)
	B = -alpha*dt - rho*sigma*(dx-mu*dt);
	C = alpha^2*dt^2 + 2*rho*sigma*alpha*dt*(dx-mu*dt) + sigma^2*(dx-mu*dt)^2 - 2*v(t)^2*a*sigma^2*(1-rho^2)*dt;
	% Equation (30) to update the variance
    if B^2 - C > 0;
        v(t+1) = sqrt(B^2 - C) - B;
    else
        % If v(t+1) is imaginary use the approximation Equation (33)
        bt = ((v(t)-alpha*dt)^2 - 2*rho*sigma*(v(t)-alpha*dt)*(dx-mu*dt) + sigma^2*(dx-mu*dt)^2)  / (2*sigma^2*(1-rho^2)*dt);
        if bt/a > 0
            % Equation (33)
            v(t+1) = sqrt(bt/a);
        else
            % If v(t+1) is still negative, take the previous value
            v(t+1) = v(t);
		end
	end
	% Equation (15) and (16)
	bt = ((v(t+1)-alpha*dt)^2 - 2*rho*sigma*(v(t+1)-alpha*dt)*(dx-mu*dt) + sigma^2*(dx-mu*dt)^2)  / (2*sigma^2*(1-rho^2)*dt);
	x1 = ((2*betap+rho*sigma*dt)*(v(t+1)-alpha*dt) - (2*rho*sigma*betap+sigma^2*dt)*(dx-mu*dt))   / (2*sigma^2*(1-rho^2)*dt);
	x2 = -2*sqrt(a*bt);
	% Combined exponent for Equation (34)
	E = exp(x1 + x2) / D;
	if method==1
		% Equation (34) for the likelihood L(t+1)
		L(t+1) = (a*bt)^(-1/4) * E * L(t);
	elseif method==2
		% Alternatively, use the log-likelihood, log of Equation (34)
		L(t+1) = -1/4*log(a*bt) + x1 + x2 -log(D) + L(t);
	end
end

% Negative likelihood is the last term.
% Since we maximize the likelihood, we minimize the negative likelihood.
y = -real(L(T));

