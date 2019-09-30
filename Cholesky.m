%Author:TAN Yueyin
%DATE:2019/9/30
function [L]=Cholesky(A)
n=size(A,1);
for k=1:n
    B(k,k)=sqrt(B(k,k));
    B(k+1:n,k)=B(k+1:n,k)/B(k,k);
    for j=k+1:n
        B(j:n,j)=B(j:n,j)-B(j:n,k)*B(j,k);
    end
end
L=tril(B);