format longEng;
clc;

k=56;
I=eye(ord);

Zdis=zeros(ord,ord);
Ydis=zeros(ord,ord);
Tidis=zeros(ord,ord);

for o=1:ord
    Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord);
    Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord);
    Tidis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
end 

gamma=g(k,2);
lamda=gamma.^2;

Q=(Ydis*Zdis)-(lamda*I)

G=Q*Tidis(:,2)

for o=1:ord
    N(o)=norm(Q(:,o));
end 

N