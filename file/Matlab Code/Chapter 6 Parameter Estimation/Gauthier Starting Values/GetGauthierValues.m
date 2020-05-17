function [sigma rho] = GetGauthierValues(kappa,theta,v0,S,K1,K2,Put1,Put2,tau,rf,q,method)

% Get the Gauthier and Possamai starting values for a single maturity
% kappa = arbitrary value for kappa
% theta = arbitrary value for theta
% v0    = arbitrary value for v0
% S  = Spot price
% K1 = First strike
% K2 = Second strike
% Put1 = Market put price for K1
% Put2 = Market put price for K2
% tau = Maturity
% rf = risk free rate
% q  = dividend yield
% Method = 1  Closed form expressions for sigma and rho
%          2  Find sigma and rho through optimization


% Find the Gauthier coefficients for each strike
[A1 B1 C1 D1] = GauthierCoefficients(kappa,theta,v0,S,K1,tau,rf,q);
[A2 B2 C2 D2] = GauthierCoefficients(kappa,theta,v0,S,K2,tau,rf,q);

if method==1
    % Closed form expressions for sigma and rho
    Rho(1) = -(-1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D2*B1+1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D1*B2-D2*A1+D2*Put1+D1*A2-D1*Put2)*2^(1/2)/((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2)/(C2*D1-C1*D2);
    Rho(2) =   (-1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D2*B1+1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D1*B2-D2*A1+D2*Put1+D1*A2-D1*Put2)*2^(1/2)/((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2)/(C2*D1-C1*D2);
    Rho(3) =  -(-1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D2*B1+1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D1*B2-D2*A1+D2*Put1+D1*A2-D1*Put2)*2^(1/2)/((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2)/(C2*D1-C1*D2);
    Rho(4) =   (-1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D2*B1+1/2*(-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2)*D1*B2-D2*A1+D2*Put1+D1*A2-D1*Put2)*2^(1/2)/((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2)/(C2*D1-C1*D2);
    Sigma(1) =  1/2*2^(1/2)*((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2);
    Sigma(2) = -1/2*2^(1/2)*((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2+(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2);
    Sigma(3) =  1/2*2^(1/2)*((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2);
    Sigma(4) = -1/2*2^(1/2)*((-C1^2*D2*B2-B1*C2^2*D1-2*D2^2*A1*B1-2*D1^2*A2*B2+2*D1^2*B2*Put2-2*D2*B1*D1*Put2+2*D2^2*B1*Put1+2*D1*A2*D2*B1+C1*C2*D2*B1+2*D1*B2*D2*A1-2*D1*B2*D2*Put1+C1*C2*D1*B2-(B1^2*C2^4*D1^2+C1^4*D2^2*B2^2+4*D2^3*A1*B1*C1^2*B2-8*D2^2*A1*B1*C1*C2*D1*B2+4*C1^2*C2^2*D2*B1*D1*B2-2*C1^3*C2*D2^2*B1*B2+4*D1^3*A2*B2*B1*C2^2-2*C1*C2^3*D2*B1^2*D1-4*D2^3*B1*Put1*C1^2*B2+4*D1*B2^2*D2^2*Put1*C1^2-4*D1^2*B2*D2*Put1*B1*C2^2-8*D1^2*B2^2*D2*Put1*C1*C2+8*D2^2*B1*Put1*C1*C2*D1*B2-4*C1^2*D2^3*A2*B1^2+4*Put1*C2^2*D1^3*B2^2-4*A1*C2^2*D1^3*B2^2+C1^2*C2^2*D1^2*B2^2+C1^2*C2^2*D2^2*B1^2+4*Put2*C1^2*D2^3*B1^2-4*B1*C2^2*D1^3*B2*Put2+4*B1^2*C2^2*D1^2*D2*Put2-4*C1^2*D2^2*B2*B1*D1*Put2+8*D1^2*B2*Put2*C1*C2*D2*B1-8*D2^2*B1^2*D1*Put2*C1*C2-4*C1^2*D2^2*B2^2*D1*A1-2*C1^3*D2*B2^2*C2*D1+8*D1*A2*D2^2*B1^2*C1*C2+4*D1*A2*D2^2*B1*C1^2*B2+8*D1^2*B2^2*D2*A1*C1*C2+4*D1^2*B2*D2*A1*B1*C2^2-4*D1^2*A2*D2*B1^2*C2^2-8*D1^2*A2*D2*B1*C1*C2*B2-2*C1*C2^3*D1^2*B2*B1)^(1/2))/(-2*D1*B2*D2*B1+D2^2*B1^2+D1^2*B2^2))^(1/2);
    % Find the allowable values for sigma and rho
    for i=1:4
        if abs(Rho(i))<1 && Sigma(i)>0;
            rho = Rho(i);
            sigma = Sigma(i);
        end
    end
elseif method==2 
    % Find the value of sigma and rho using constrained optimization
    e = 1e-5;
    lb = [e  -.999];  % Lower bound on the estimates
    ub = [10  .999];  % Upper bound on the estimates
    start = [.3 -0.9];
    Coeff1 = [A1 B1 C1 D1];
    Coeff2 = [A2 B2 C2 D2];
    [SigmaRho] = fmincon(@(p) GauthierObjFun(p,Coeff1,Coeff2,Put1,Put2),start,[],[],[],[],lb,ub);
    sigma = SigmaRho(1);
    rho   = SigmaRho(2);
end

