function [Ti,g]=simple_QR_decomp(ord,num_files,Z,Y)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % ������� ��� ������ Z' ��� �� ������������ ��������� - (ord x ord)
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % ������� ��� ������ �' ��� �� ������������ ��������� - (ord x ord)
    end
    
    [Tidis,ldis]=eig(Ydis*Zdis); % ����������� ��������������� ��� ��������� �� ���� ��� QR decomposition - ��� �� ��� ������� (ord x ord)
    gdis=sqrt(diag(ldis)); % ����������� ��� propagation constant - (1 x ord)
    
    g(k,:)=gdis; % ���������� ��� propagation constant - (num_files x ord)
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % ���������� ��� ������ ��������������� - (num_files x ord^2)
    end
    
end