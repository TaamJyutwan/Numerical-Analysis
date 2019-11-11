%Author: TAN Yueyin
%DATE: 2019/11/10
function [v]=EigenValue_OriginShift(A)
n = size(A,1);
A = Hessenberg(A);
i = 0;
while 1
    for t = 2:n
        if abs(A(t-1,t-1) * A(t,t) -det(A(t-1:t,t-1:t))) < 1e-10
            A(t,t-1)=0;
        end
    end
    %将次对角线上对特征值影响较小的值变为0
    for i = n:-1:2
        if i ~= 2 && A(i,i-1) ~= 0 && A(i-1,i-2) ~= 0
            miu = A(n,n);
            A1 = A - miu * eye(n);
            [Q,R] = qr(A1);
            A = R * Q + miu * eye(n);
            break
        end
    end
    if i == 2 %若上一步循环能够到i=2，说明全部复特征值对应的2*2子块已经显现
        break
    end  
end
for i = 2:n
    if A(i,i-1) ~= 0
        p = A(i,i) + A(i-1,i-1);
        q = A(i,i) * A(i-1,i-1) - A(i,i-1) * A(i-1,i);
        delta = p ^ 2 - 4 * q;
        A(i-1,i-1) = (p + sqrt(delta)) / 2;
        A(i,i) = (p-sqrt(delta)) / 2;
    end
end
v = diag(A);