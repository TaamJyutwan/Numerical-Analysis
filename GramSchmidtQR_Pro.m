function [Q,R]=GramSchmidtQR_Pro(A)
m=size(A,1);
n=size(A,2);
Q=zeros(m,n);
R=zeros(n);
for j=1:n
    y=A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*y;
        y=y-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(y,2);
    Q(:,j)=y/R(j,j);
end