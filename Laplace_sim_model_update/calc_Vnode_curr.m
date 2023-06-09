function [V]=calc_Vnode_curr(ord,f,Ybranch,Io1)

V=zeros(2*ord,max(size(f))); % (2*ord x max(size(f)))
Ybranch_dis_f=zeros(2*ord,2*ord); % (2*ord x 2*ord)
for k=1:1:max(size(f))
    
    for o=1:1:2*ord
        Ybranch_dis_f(o,:)=Ybranch(k,(o-1)*2*ord+1:o*2*ord); % ��������� ��� ������ Ybranch ��� ������������� ���������� ��� �������� �� ����������� (2*ord x 2*ord) ���� �� �������������� �� ������� �������
    end
    
    Ired=[Io1(k); zeros(length(Ybranch_dis_f)-1,1)];
    Vred=Ybranch_dis_f\Ired;

    V(:,k)=Vred; % ���������� ���� ��� ������ ��� ������ �� ��� �� ����� ���������� - (2*ord x max(size(f)))
end