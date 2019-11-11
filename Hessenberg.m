function [A]=Hessenberg(A)
n = size(A,1);
for k = 1 : n - 2
    [v,b] = Householder(A(k+1:n,k));
    h = eye(n-k) - b * (v * v');
    A(k+1:n,k:n) = h * A(k+1:n,k:n);
    A(1:n,k+1:n) = A(1:n,k+1:n) * h;
end

    
    
