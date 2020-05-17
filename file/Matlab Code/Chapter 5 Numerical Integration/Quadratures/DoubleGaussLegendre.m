function y = DoubleGaussLegendre(f,a,b,c,d,x1,w1,x2,w2)

% Double integration by composite Gauss Legendre rule
% f = bivariate function to be integrated
% x = vector of GL abscissas
% w = vector of GL weights
% Int(x:a,b; y:c,d);


N1 = length(x1);
N2 = length(x2);
h1 = (b-a)/2;
h2 = (b+a)/2;
k1 = (d-c)/2;
k2 = (d+c)/2;

for i=1:N1;
    for j=1:N2
        Int(i,j) = h1*k1*w1(i)*w2(j) * f(h1*x1(i)+h2,k1*x2(j)+k2);
    end
end
y = sum(sum(Int));
