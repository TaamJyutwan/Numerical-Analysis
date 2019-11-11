function [H] = FrancisQR(H)
n = size(H,1);
m = n - 1;
s = H(m,m) + H(n,n);
t = H(m,m) * H(n,n) - H(m,n) * H(n,m);
x = H(1,1) * H(1,1) + H(1,2) * H(2,1) - s * H(1,1) + t;
y = H(2,1) * (H(1,1) + H(2,2) - s);
z = H(2,1) * H (3,2);
for k = 0 : n - 3
    [v,b] = Householder([x,y,z]');
    h = eye(3) - b * (v * v');
    q = max(1,k);
    H(k+1:k+3,q:n) = h * H(k+1:k+3,q:n);
    r = min(k+4,n);
    H(1:r,k+1:k+3) = H(1:r,k+1:k+3) * h;
    x = H(k+2,k+1);
    y = H(k+3,k+1);
    if k < n-3
        z = H(k+4,k+1);
    end
end
[v,b] = Householder([x,y]');
h = eye(2) - b * (v * v');
H(n-1:n,n-2:n) = h * H(n-1:n,n-2:n);
H(1:n,n-1:n) = H(1:n,n-1:n) * h;

