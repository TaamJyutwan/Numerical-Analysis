%Author: TAN Yueyin
%DATE: 2019/10/29
%用于解 A * x = b 的实用共轭梯度法
%方程的系数矩阵 A 为对称正定矩阵
function [x]=PragmaticConjugateGradient(A,b,km,acr)
n = size(A,1);
x = rand(n,1);
r = b - A * x;
k = 0;
l = r' * r;
while norm(r,2) > acr && k < km
    k=k+1;
    if k == 1
        p = r;
    else
        beta = l / lp;
        p = r + beta * p;
    end
    w = A * p;
    a = l / (p' * w);
    x = x + a * p;
    r = r - a * w;
    lp = l;
    l = r' * r;
end