function [x] = GaussSeidel_Iteration(A,b,acr)
n = size(A,1);
x = rand(n,1);
while 1
    xs = x;
    for i = 1 : n
        x(i) = ( -A(i,1:i-1) * x(1:i-1) - A(i,i+1:n) * x(i+1:n) + b(i)) / ( A(i,i));
    end
    if norm(x - xs,2) < acr
        break
    end
end
    
