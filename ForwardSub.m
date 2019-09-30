%Author:TAN Yueyin
%DATE:2019/9/30
function [b]=ForwardSub(L,b)
n=size(L,1);
for i=1:n-1
    b(i)=b(i)/L(i,i);
    b(i+1:n)=b(i+1:n)-b(i)*L(i+1:n,i);
end
b(n)=b(n)/L(n,n);