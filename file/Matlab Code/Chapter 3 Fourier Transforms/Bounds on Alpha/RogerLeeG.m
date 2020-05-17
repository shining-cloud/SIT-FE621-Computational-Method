function y = RogerLeeG(phi,kappa,rho,sigma)

% function to find g for Roger Lee bounds

A = (rho*sigma*phi*i - kappa);
B = sigma^2*(phi*i + phi^2);
d = sqrt(A^2 + B);
g = (kappa - rho*sigma*phi*i + d) / (kappa - rho*sigma*phi*i - d);

y = g;


