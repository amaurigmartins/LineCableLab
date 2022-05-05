function [Ti,g]=SQP_calc_norm_back(ord,num_files,Z,Y,freq)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

O=eye(ord); % Unit matrix - (ord x ord)

e0=8.85418782*1e-12; % Permittivity �����-����
m0=4*pi*1e-7; % Permeability �����-����

for k=num_files:-1:1 % �� ������ ��� ����������!
    k % ������� ��� �� ���������� ��� command window �� ���� ��������� �����������
    
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % ������� ��� Z' ��� �� ������� ��������� - (ord x ord)
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % ������� ��� �' ��� �� ������� ��������� - (ord x ord)
    end 
        
    if k==num_files || k==1 % ���� ����� ��� ��������� ��������� ����������� QR decomposition! ���� ����� ����� ���� ������ � ���������� ��� ���� ��������� ����� ����� �� 0 Hz ��� ��� ������ �� ��������������� �� ����� �� ��� ������� �������� - �������!!! ����������� ������ �� ���������� ��� ��� ������� ��������� �.�.�.
        S=Ydis*Zdis; % �������� �'�' - (ord x ord)
        [Tidis,ldis]=eig(S); % QR decomposition - ���� �� ������� ����� (ord x ord)
        gdis=sqrt(diag(ldis)); % ������ ��� propagation constant - (1 x ord)

        g(k,:)=gdis; % ���������� ��� propagation constant - (num_files x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % ���������� ��� ������ ��������������� - (num_files x ord^2)
        end

    else % ���� ��������� ���������� ����������� � ���������� SQP
        S=Ydis*Zdis; % �������� �'�' - (ord x ord)
        norm=-((2*pi*freq(k)).^2)*e0*m0; % ����� � ����� ������������� �� �������� �'�' ���� �� ������������� ���������� �� ��������� - �������!!! �� ��������� ��� �� ������������ ����������� ���������� �������� �� ��� ���� ����� ��� �� �������������� ���������� ���������� ���� ����� ��������������� �� ��� ������� ������
        %S=S./norm; % �������������� ��� ��������� �'�' - (ord x ord) 
        S=(S./norm)-O;
        
        gdis=g(k+1,:); % �������� ��� propagation constant ��� ������������ ���������� - Seeding!! - (1 x ord)
        ldis=gdis.^2; % ����������� ��� ��������� - (1 x ord)
        %ldis=ldis./norm; % �������������� ��� ���������
        ldis=(ldis./norm)-1;
        
        for o=1:ord
            Tidis(o,:)=Ti(k+1,(o-1)*ord+1:o*ord); % �������� ��� ������ ��������������� ��� ������������ ���������� - Seeding!! - (ord x ord)
        end
                
        for o=1:ord
            x0=[real(Tidis(:,o));imag(Tidis(:,o));real(ldis(o));imag(ldis(o))]; % �� �������� x �������� ����� �� ���������� ��� ���� �� ���������� ����� ���� ���������������, ��� ������ �� ���������� ��� ���������� ����� ��� ���������� ��������� - (2*ord+2 x 1)
            Q=[real(S) -imag(S);imag(S) real(S)]; % ���������� ������� ��� ��� ��������� ��� ���������� ������������ (2*ord x 2*ord)
            obj=@(x)((transpose(Q*x(1:2*ord)-[x(2*ord+1)*O -x(2*ord+2)*O;x(2*ord+2)*O x(2*ord+1)*O]*x(1:2*ord)))*(Q*x(1:2*ord)-[x(2*ord+1)*O -x(2*ord+2)*O;x(2*ord+2)*O x(2*ord+1)*O]*x(1:2*ord)))/2; % ��������� ������������ - ����� ������� (scalar)! ���������� ���� ���������� ������
            
            %c=@(x)[(2*pi*freq(k))/imag(sqrt(((x(2*ord+1)+1i*x(2*ord+2))+1)*norm));-(2*pi*freq(k))/imag(sqrt(((x(2*ord+1)+1i*x(2*ord+2)+1)*norm)));-real((sqrt((x(2*ord+1)+1i*x(2*ord+2)+1))*norm))]-[1/sqrt(e0*m0);0;0]; % �� ��������� ����������� ��� ������ c<=0 - (3 x 1)
            c=@(x)[];
            ceq=@(x)[sum((x(1:ord).^2)-(x(ord+1:2*ord).^2));sum(x(1:ord).*x(ord+1:2*ord))]-[1;0]; % �� ��������� ����������� ��� ������ ceq=0 - (2 x 1)
            
            nonlincon=@(x)deal(c(x),ceq(x)); % ���������� ��� �� ��������� ����������� �� ��� ��������� �� �������� ������ �� �������� x
            options=optimset('Algorithm','sqp','Display','off','FinDiffType','central','MaxFunEvals',6000); % �������� ���������� ��� ��� ������� ��� ����������� ��� ��������������� ��� ����������-������������
            
            x=fmincon(obj,x0,[],[],[],[],[],[],nonlincon,options); % �������������� ��� ���������� ������������ - ��������� ��� ����������� x ��� ����������� ���� ��� �������������� - (2*ord+2 x 1)
            
            Tidis(:,o)=x(1:ord)+1i*x(ord+1:2*ord); % ������� ��� ������ ��������������� - (ord x ord)
            ldis(o)=x(2*ord+1)+1i*x(2*ord+2); % ������� ��� ����������� ��������� - (ord x ord)
        end
        
        %ldis=ldis.*norm; % ���������� ��� ��������������� - (1 x ord)
        ldis=(ldis+1).*norm;
        g(k,:)=sqrt(ldis); % ����������� ��� ���������� ��� propagation constant - (1 x ord)
                
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % ���������� ��� ������ ��������������� - (num_files x ord^2)
        end
        
    end
end
            
        