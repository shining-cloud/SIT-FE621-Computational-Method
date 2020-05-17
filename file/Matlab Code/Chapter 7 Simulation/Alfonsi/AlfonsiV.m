function newV = AlfonsiV(param,vt,dt);

kappa = param(1);
theta = param(2);
sigma = param(3);
v0    = param(4);
rho   = param(5);

phi = (1-exp(-kappa*dt/2))/kappa;
S = (sigma^2/4 - theta*kappa);
E = exp(kappa*dt/2);

% K2 parameter
if sigma^2 > 4*kappa*theta
	K2 = E*(S*phi + (sqrt(E*S*phi) + sigma/2*sqrt(3*dt))^2);
else
	K2 = 0;
end

% Update the variance
if vt >= K2
	U = rand(1);
	if U <= 1/6;
		Y = sqrt(3);
	elseif U <= 1/3;
		Y = -sqrt(3);
	else
		Y = 0;
    end
	phi = (1-exp(-kappa*dt/2))/kappa;
	S = (theta*kappa - sigma^2/4);
	E = exp(-kappa*dt/2);
	newV = E*(sqrt(S*phi + E*vt) + sigma/2*sqrt(dt)*Y)^2 + S*phi;
else
	[u1 u2] = CIRmoments(param,vt,dt);
	Pi = 0.5 - 0.5*sqrt(1 - u1^2/u2);
	U = rand(1);
	if U <= Pi
		newV = u1/2/Pi;
	elseif U > Pi
		newV = u1/2/(1-Pi);
    end
end

