function [Ybranch_dis]=get_Ybranch(ord,length,Ti_vector,Ys,Yr,gamma,Z_vector)

Z=zeros(ord,ord);
Ti=zeros(ord,ord);

for o=1:1:ord
    Z(o,:)=Z_vector(1,(o-1)*ord+1:o*ord); % ��������� ��� ������ Z' ��� ������������� ��������� ���������� �� ����������� (ord x ord) ���� �� �������������� �� ������� �������
    Ti(o,:)=Ti_vector(1,(o-1)*ord+1:o*ord); % ��������� ��� ������ Ti ��� ������������� ��������� ���������� �� ����������� (ord x ord) ���� �� �������������� �� ������� �������
end

Tv=transpose(inv(Ti)); % ����������� ��� ����������� ������ Tv ��� �� ������������ �������� ��������� - (ord x ord)

%B=Tv*diag(sinh(gamma*length))*inv(diag(gamma))*inv(Tv)*Z; % ����������� ��� ������ B ��� �� ������������ �������� ��������� - (ord x ord)
B=Tv*(diag(sinh(gamma*length))/(diag(gamma)))*(Tv\Z);
%B=Tv*(diag((exp(gamma*length)-exp(-gamma*length))/2)/(diag(gamma)))*(Tv\Z);

%D=Tv*diag(cosh(gamma*length))*inv(Tv); % ����������� ��� ������ D ��� �� ������������ �������� ��������� - (ord x ord)
D=Tv*(diag(cosh(gamma*length))/Tv);
%D=Tv*(diag((exp(gamma*length)+exp(-gamma*length))/2)/Tv);

%Yseries=inv(B); % ����������� ��� ������ Yseries ��� �� ������������ �������� ��������� - (ord x ord)
Yseries=B\(eye(ord)); %TEST
%Yshunt=inv(B)*D-inv(B); % ����������� ��� ������ Yshunt ��� �� ������������ �������� ��������� - (ord x ord)
%Yshunt=(B\D)-inv(B);
Yshunt=(B\D)-Yseries; %TEST


Yself=Yseries+Yshunt; % ����������� ��� ������ Yself ��� �� ������������ �������� ��������� - (ord x ord)
Ymutual=-Yseries; % ����������� ��� ������ Ymutual ��� �� ������������ �������� ��������� - (ord x ord)

Yline=[Yself Ymutual;Ymutual Yself]; % ����������� ��� ������ Yline ��� �� ������������ �������� ��������� - ((2*ord) x (2*ord))

Yload=[diag(Ys) zeros(ord);zeros(ord) diag(Yr)]; % ����������� ��� ������ Yload ��� �� ������������ �������� ��������� - ((2*ord) x (2*ord))

Ybranch=Yline+Yload; % ����������� ��� ������ Ybranch ��� �� ������������ �������� ��������� - ((2*ord) x (2*ord))

for o=1:1:2*ord
    Ybranch_dis(1,(o-1)*2*ord+1:o*2*ord)=Ybranch(o,:); % ���������� ��� ������ Ybranch ��� �� ������������ �������� ��������� - (1 x (2*ord)^2) 
end