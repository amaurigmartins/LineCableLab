function [w, Z, Y] = get_ZY_from_cyz(filename)
D=importdata(filename);     C=D(2);           N=D(3);
P=C*C*2*2;         last=length(D);      D=D(4:last,1);
last=length(D);   D=D((P+2):last);
N=N-2;     zi = sqrt(-1);     n=1;
w=zeros(N,1);   Z=zeros(C,C,N);   Y=zeros(C,C,N);
for f=1:N
    w(f)=2*pi*D(n);     q=1;
    for k=1:C
        for m=1:C
            Z(k,m,f) = (D(n+q)+ zi*D(n+q+1));
            Y(k,m,f) = (D(n+q+2)+zi*D(n+q+3));
            q = q+4;
        end
    end
    n=n+P+1;
end

end