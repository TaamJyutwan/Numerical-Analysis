%Author: TAN Yueyin
%DATE: 2019/10/29
%用于解 A * x = b 的最速下降法
%方程的系数矩阵 A 为对称正定矩阵
function [x] = SteepestDecent(A,b,acr)
n = size(A,1);
x = rand(n,1);
r = b - A * x;
k = 0;
while norm(r,2) > acr
    k = k + 1;
    a = r' * r / (r' * A * r);
    x = x + a * r;
    r = b - A * x;
end
