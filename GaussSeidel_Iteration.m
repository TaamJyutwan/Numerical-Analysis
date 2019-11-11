function [x] = GaussSeidel_Iteration(A,b,acr)
DL = tril(A);
U = - triu(A,1);
B = DL \ U;
g = DL \ b;
xs = b;
x = B * xs + g;
while norm(x - xs,2) > acr
    xs = x;
    x = B * x + g;
end