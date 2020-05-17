function y = GauthierObjFun(param,Coeff1,Coeff2,Put1,Put2)

% Returns the objective function for calculating the Gauthier and Rivaille
% coefficients under optimization

sigma = param(1);
rho   = param(2);

A1 = Coeff1(1);
B1 = Coeff1(2);
C1 = Coeff1(3);
D1 = Coeff1(4);

A2 = Coeff2(1);
B2 = Coeff2(2);
C2 = Coeff2(3);
D2 = Coeff2(4);

y = (A1 + B1*sigma^2 + C1*rho*sigma + D1*rho^2*sigma^2 - Put1)^2 ...
  + (A2 + B2*sigma^2 + C2*rho*sigma + D2*rho^2*sigma^2 - Put2)^2;

	
