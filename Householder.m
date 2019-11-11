%Author:TAN Yueyin
%DATE: 2019/11/11
function [v,b]=Householder(x)
n=length(x);
v=zeros(n,1);
x=x/norm(x,inf);
v(2:n)=x(2:n);
M=x(2:n)'*x(2:n);
if M==0
    b=0;
else
    a=sqrt(x(1)*x(1)+M);
    if x(1)<=0
        v(1)=x(1)-a;
    else
        v(1)=-M/(x(1)+a);
    end
    b=2*v(1)*v(1)/(M+v(1)*v(1));
    v=v/v(1);
end

