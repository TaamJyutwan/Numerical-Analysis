%Author: TAN Yueyin
%DATE:2019/11/10
function [v] = EigenValue_DoubleShift(A)
n = size(A,1);
H = Hessenberg(A);

while 1
    i1 = 0;
    for i = 2:n
        if abs(H(i-1,i-1) * H(i,i) -det(H(i-1:i,i-1:i))) < 1e-10
            H(i,i-1)=0;
        end
    end
    
    for i = n:-1:2
            if i ~=2 && H(i,i-1) ~= 0 && H(i-1,i-2) ~= 0
               i2 = i;
               for j = i2:-1:2
                   if H(j,j-1) == 0 %找到分块位置，可以这个元为界限划四分块
                       i1 = j;
                       H(i1:i2,i1:i2) = FrancisQR(H(i1:i2,i1:i2));
                       break
                   end
               end
               if i1 == 0
                   H(1:i2,1:i2) = FrancisQR(H(1:i2,1:i2));
               end
               break
            end
    end
    if i == 2
        break
    end
end

for i = 2:n
    if H(i,i-1) ~= 0
        p = H(i,i) + H(i-1,i-1);
        q = H(i,i) * H(i-1,i-1) - H(i,i-1) * H(i-1,i);
        delta = p ^ 2 - 4 * q;
        H(i-1,i-1) = (p + sqrt(delta)) / 2;
        H(i,i) = (p-sqrt(delta)) / 2;
    end
end
v = diag(H);