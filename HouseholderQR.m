function [Q,R]=HouseholderQR(A)
m=size(A,1);
n=size(A,2);
Q=eye(m);
for j=1:n
    if j<m
        [v,b]=Householder(A(j:m,j));
        H=eye(m-j+1)-b*(v*v');
        h=[eye(j-1), zeros(j-1,m-j+1);zeros(m-j+1,j-1),H];
        A=h*A;
        Q=Q*h;
    end
end
R=triu(A);
