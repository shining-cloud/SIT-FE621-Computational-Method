function y = DoubleTrapz(f,X,Y)

% Double integration by Composite Trapezoidal rule
% f = bivariate function to be integrated
% X = vector for the X-direction
% Y = vector for the Y-direction

Nx = length(X);
Ny = length(Y);

for y=2:Ny
    a = Y(y-1);
    b = Y(y);
    for x=2:Nx
        c = X(x-1);
        d = X(x);
        term1 = f(a,c) + f(a,d) + f(b,c) + f(b,d);
        term2 = f((a+b)/2,c) + f((a+b)/2,d) + f(a,(c+d)/2) + f(b,(c+d)/2);
        term3 = f((a+b)/2,(c+d)/2);
        Int(x,y) = (b-a)*(d-c)/16*(term1 + 2*term2 + 4*term3);
    end
end

y = sum(sum(Int));
