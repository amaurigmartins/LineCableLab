function [Ti,g]=product_correl_back(ord,num_files,Z,Y,freq)

[Ti,g]=simple_QR_decomp(ord,num_files,Z,Y); % ���� QR decomposition ����� �������� ��� eigenvector switchover - ���������� ������ ��������������� Ti (current) (num_files x ord^2) ��� ��������� g (sqrt) (num_files x ord)

%[Ti,g]=atp_Tsiam(Z,Y,ord,freq,num_files); % � ��������� ��� ���������� - ����������� ��� QR decomposition ��� �������������� ��� ����������� ������ �� ���� ��� ����������� ��� ATP/EMTP - ���������� ������ ��������������� Ti (current) ((num_files x ord^2)) ��� ��������� g (sqrt) (num_files x ord)

Ti_prev=zeros(ord,ord); % O ������� ��������������� ��� �������� ���������� - ����������� ��� switchover effect! - (ord x ord)
Ti_next=zeros(ord,ord); % O ������� ��������������� ��� �������� ���������� - �� ����������� ��� switchover effect! - (ord x ord)
flag=0; % ������� ��� ��� ����� ���������� ��� ������ ���������������

for k=num_files:-1:2 % �� ������ ��� ����������!
        
    for o=1:ord
        Ti_prev(o,:)=Ti(k,(o-1)*ord+1:o*ord); % ������� ��� ������ ��������������� ��� �������� ���������� - ����������� ��� switchover effect! - (ord x ord)
        Ti_next(o,:)=Ti(k-1,(o-1)*ord+1:o*ord); % ������� ��� ������ ��������������� ��� �������� ���������� - �� ����������� ��� switchover effect! - (ord x ord)
    end
    
    sign=1; % ������� ��� �� loop ��� ��������� ��� eigenvector switchover
    while sign==1
        sign=0; % Reset ��� �������� - �� ��� ������� ������ ������ ���� �� loop �� ����������
        Ti_next_mod=ctranspose(Ti_next); % ����������� ��� complex conjugate transpose ������ ��������������� ��� �������� ���������� - (ord x ord)
        P=Ti_next_mod*Ti_prev; % ����������� ��� ������������� ��������� - (ord x ord)
    
        for o=1:ord % ��� ���� ������ ��� ������������� ������ ���������
            [~,i]=max(P(o,:)); % ����������� ��� ������ ���� ��������� �� ������� �������� ��� ������� "�" ��� ������������� ������ ���������
            if i~=o % ������ � ������� ��� ������ ����� �������� ��� �� ������ ��� �������
                sign=1; % ������� ��� �� ���������� ��� loop ��� ��������� ���� ������� ������ ��� ��������
                flag=1; % ������� ��� ��������� ���������� ���� ������� ������ ��� ��������
                
                temp=Ti_next(:,i);
                Ti_next(:,i)=Ti_next(:,o); % �������� ��� switchover ���� ������ ���������������
                Ti_next(:,o)=temp;
            
                temp=g(k-1,i);
                g(k-1,i)=g(k-1,o); % �������� ��� switchover ���� ������ ��������� ��� ����� ���������� - (num_files x ord)
                g(k-1,o)=temp;
                break % ������ ��� �� for ���� �� ��������������� �� ������������ �������� �� ��� �� ����� ���������� ������ ���������������� ��� �������� ����������
            end
        end
    end
    
    if flag==1
        flag=0; % Reset ��� ��������
        k-1 % ������� ��� command window ��� ������ switchover �� ���� �� ���������
        
        for o=1:ord
            Ti(k-1,(o-1)*ord+1:o*ord)=Ti_next(o,:); % ����� ���������� ��� ������ ��������������� ���� ������������ ����� ��������� ���� �� ������ ��� - (num_files x ord^2)
        end
    end
end

