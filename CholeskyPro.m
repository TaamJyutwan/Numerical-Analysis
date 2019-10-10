%Author:TAN Yueyin
%DATE:2019/9/30
function [L,D]=CholeskyPro(A)
n=size(A,1);
v=zeros(n,1);
for i=1:n
    for j=1:i-1
        v(j,1)=A(i,j)*A(j,j);
    end
    A(i,i)=A(i,i)-A(i,1:i-1)*v(1:i-1,1);
    A(i+1:n,i)=(A(i+1:n,i)-A(i+1:n,1:i-1)*v(1:i-1,1))/A(i,i);
end
L=tril(A,-1)+eye(n);
D=diag(diag(A));