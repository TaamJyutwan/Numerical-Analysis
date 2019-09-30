%Author:TAN Yueyin
%DATE:2019/9/30
function [b]=BackwardSub(U,b)
n=size(U,1);
for i=n:-1:2
    b(i)=b(i)/U(i,i);
    b(i-1:-1:1)=b(i-1:-1:1)-b(i)*U(i-1:-1:1,i);
end
b(1)=b(1)/U(1,1);