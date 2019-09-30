%Author:TAN Yueyin
%DATE:2019/9/30
function [L,D]=CholeskyPro(A)
n=size(A,1);
for i=1:n
    for j=1:i-1
        v(j,1)=B(i,j)*B(j,j);
    end
    B(i,i)=B(i,i)-B(i,1:i-1)*v(1:i-1,1);
    B(i+1:n,i)=(B(i+1:n,i)-B(i+1:n,1:i-1)*v(1:i-1,1))/B(i,i);
end
L=tril(B,-1)+eye(n);
D=diag(diag(B));