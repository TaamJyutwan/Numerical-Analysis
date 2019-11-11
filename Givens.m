function [c,s]=Givens(a,b)
if b==0
    c=1;
    s=0;
else
    if abs(b)>abs(a)
        t=a/b;
        s=1/sqrt(1+t*t);
        c=t*s;
    else
        t=b/a;
        c=a/sqrt(1+t*t);
        s=t*c;
    end
end
