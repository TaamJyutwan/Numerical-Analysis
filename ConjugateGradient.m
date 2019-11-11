%Author: TAN Yueyin
%DATE: 2019/10/29
%用于解 A * x = b 的共轭梯度法
%方程的系数矩阵 A 为对称正定矩阵
function [x]=ConjugateGradient(A,b,acr)
n = size(A,1);
x = rand(n,1);
r = b - A * x;
k = 0;
while norm(r,2) > acr
    k = k + 1;
    if k == 1
        p = r;
    else
        beta = r' * r /(rs' * rs);
        p = r + beta * p;
    end
    a = r' * r / (p' * A * p);
    x = x + a * p;
    rs = r;
    r = r - a * A * p;
end
   