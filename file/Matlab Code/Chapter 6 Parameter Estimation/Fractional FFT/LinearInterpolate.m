function Yi = LinearInterpolate(X,Y,Xi)

N = length(X);
M = length(Xi);
for j=1:M
    k = find(Xi(j)<=X);
    k = k(1)-1;
    Yi(j) = Y(k+1)*(Xi(j)-X(k))/(X(k+1)-X(k)) + Y(k)*(X(k+1)-Xi(j))/(X(k+1)-X(k));
end


