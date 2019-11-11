%Author: TAN Yueyin
%DATE: 2019/11/10
function [x] = Jacobi_Iteration(A,b,acr)
D = diag(diag(A));
L = - tril(A,-1);
U = - triu(A,1);
B = D \ (L+U);
g = D \ b;
x = b;
while norm(x-B*x-g) > acr
    x = B * x + g;
end

