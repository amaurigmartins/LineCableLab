function [Ybranch_dis]=get_Ybranch_scaled(ord,length,Ti_vector,Ys,Yr,gamma,Z_vector)

Z=zeros(ord,ord);
Ti=zeros(ord,ord);

for o=1:1:ord
    Z(o,:)=Z_vector(1,(o-1)*ord+1:o*ord); % ��������� ��� ������ Z' ��� ������������� ��������� ���������� �� ����������� (ord x ord) ���� �� �������������� �� ������� �������
    Ti(o,:)=Ti_vector(1,(o-1)*ord+1:o*ord); % ��������� ��� ������ Ti ��� ������������� ��������� ���������� �� ����������� (ord x ord) ���� �� �������������� �� ������� �������
end

Tv=transpose(inv(Ti)); % ����������� ��� ����������� ������ Tv ��� �� ������������ �������� ��������� - (ord x ord)

h=exp(-gamma*length);

alpha=(gamma.*(1+h.^2))./(1-h.^2);
beta=(-2*gamma.*h)./(1-h.^2);

%A=inv(Z)*Tv*diag(alpha)*inv(Tv);
A=(Z\eye(ord))*Tv*diag(alpha)*(Tv\eye(ord));

%B=inv(Z)*Tv*diag(beta)*inv(Tv);
B=(Z\eye(ord))*Tv*diag(beta)*(Tv\eye(ord));

Yline=[A B;B A]; % ����������� ��� ������ Yline ��� �� ������������ �������� ��������� - ((2*ord) x (2*ord))

Yload=[diag(Ys) zeros(ord);zeros(ord) diag(Yr)]; % ����������� ��� ������ Yload ��� �� ������������ �������� ��������� - ((2*ord) x (2*ord))

Ybranch=Yline+Yload; % ����������� ��� ������ Ybranch ��� �� ������������ �������� ��������� - ((2*ord) x (2*ord))

for o=1:1:2*ord
    Ybranch_dis(1,(o-1)*2*ord+1:o*2*ord)=Ybranch(o,:); % ���������� ��� ������ Ybranch ��� �� ������������ �������� ��������� - (1 x (2*ord)^2) 
end