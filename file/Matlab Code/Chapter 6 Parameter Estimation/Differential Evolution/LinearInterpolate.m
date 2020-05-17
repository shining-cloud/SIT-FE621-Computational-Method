function Yi = LinearInterpolate(X,Y,Xi)

% Function for linear documentation
% X = vector of X's
% Y = vector of Y's
% Xi = value of X where Y should be interpolated

N = length(X);
M = length(Xi);
for j=1:M
    k = find(Xi(j)<=X);
    k = k(1)-1;
    Yi(j) = Y(k+1)*(Xi(j)-X(k))/(X(k+1)-X(k)) + Y(k)*(X(k+1)-Xi(j))/(X(k+1)-X(k));
end


