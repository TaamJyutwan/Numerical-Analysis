%Author:TAN Yueyin
%DATE:2019/9/30
function [L,U,b]=GaussColumn(A,b)
n=size(A,1);
Ab=[A,b];
for k=1:n-1
    MaxElement=max(abs(A(k:n,k)));
    [p,q]=find(abs(A)==MaxElement);
    Ab([k,p],:)=Ab([p,k],:);
    if A(k,k)==0
        stop
    else
        Ab(k+1:n,k)=Ab(k+1:n,k)/Ab(k,k);
        Ab(k+1:n,k+1:n)=Ab(k+1:n,k+1:n)-Ab(k+1:n,k)*Ab(k,k+1:n);
    end
end
L=tril(Ab(1:n,1:n),-1)+eye(n);
U=triu(Ab(1:n,1:n));
b=Ab(1:n,n+1);
    