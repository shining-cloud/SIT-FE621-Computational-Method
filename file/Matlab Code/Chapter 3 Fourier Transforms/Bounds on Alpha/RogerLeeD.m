function y = RogerLeeD(phi,kappa,rho,sigma)

% Function to find D for Roger Lee bounds

A = (rho*sigma*phi*i - kappa);
B = sigma^2*(phi*i + phi^2);
d = sqrt(A^2 + B);

y = d;


