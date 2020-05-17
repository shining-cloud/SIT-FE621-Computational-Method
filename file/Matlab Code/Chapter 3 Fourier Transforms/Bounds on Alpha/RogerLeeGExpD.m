function y = RogerLeeGExpD(phi,kappa,rho,sigma,tau)

% Function to find g*exp(d) for Roger Lee bounds

A = (rho*sigma*phi*i - kappa);
B = sigma^2*(phi*i + phi^2);
d = sqrt(A^2 + B);
g = (kappa - rho*sigma*phi*i + d) / (kappa - rho*sigma*phi*i - d);

E = real(g*exp(d*tau));

y = (E - 1)^2;



