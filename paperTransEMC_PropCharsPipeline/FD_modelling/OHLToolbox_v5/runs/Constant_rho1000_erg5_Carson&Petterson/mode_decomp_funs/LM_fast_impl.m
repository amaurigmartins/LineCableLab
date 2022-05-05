function [Ti,g]=LM_fast_impl(ord,num_files,Z,Y,freq)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

O=eye(ord); % Unit matrix - (ord x ord)
error=1e-8; % ���� ��������� ���� ���������� � ����������

e0=8.85418782*1e-12; % Permittivity �����-����
m0=4*pi*1e-7; % Permeability �����-����

for k=1:1:num_files % �� ������ ��� ����������!
    k
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % ������� ��� Z' ��� �� ������� ��������� - (ord x ord)
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % ������� ��� Y' ��� �� ������� ��������� - (ord x ord)
    end   
    
    if k==1 % ���� ����� ��������� ����������� QR decomposition! - �������!!! ����������� ������ �� ���������� ��� ��� ������� ��������� �.�.�.
        S=Ydis*Zdis; % �������� Y'Z' - (ord x ord)
        [Tidis,ldis]=eig(S); % QR decomposition - ���� �� ������� ����� (ord x ord)
        gdis=sqrt(diag(ldis)); % ������ ��� propagation constant - (1 x ord)

        g(k,:)=gdis; % ���������� ��� propagation constant - (num_files x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % ���������� ��� ������ ��������������� - (num_files x ord^2)
        end

    else % ���� ��������� ���������� ����������� � ���������� Newton-Raphson
        S=Ydis*Zdis; % �������� Y'Z' - (ord x ord)
        norm=-((2*pi*freq(k)).^2)*e0*m0; % ����� � ����� ������������� �� �������� Y'Z' ���� �� ������������� ���������� �� ��������� - �������!!! �� ��������� ��� �� ������������ ����������� ���������� �������� �� ��� ���� ����� ��� �� �������������� ���������� ���������� ���� ����� ��������������� �� ��� ������� ������
        %S=S./norm; % �������������� ��� ��������� Y'Z' - (ord x ord) 
        S=(S./norm)-O;
        
        gdis=g(k-1,:); % �������� ��� propagation constant ��� ������������ ���������� - Seeding!! - (1 x ord)
        ldis=gdis.^2; % ����������� ��� ��������� - (1 x ord)
        %ldis=ldis./norm; % �������������� ��� ���������
        ldis=(ldis./norm)-1;

        for o=1:ord
            Tidis(o,:)=Ti(k-1,(o-1)*ord+1:o*ord); % �������� ��� ������ ��������������� ��� ������������ ���������� - Seeding!! - (ord x ord)
        end
        
        for o=1:ord
            x_old=[Tidis(:,o);ldis(o)]; % �� �������� x �������� ��� ������������ ��� ��� ���������� �������� - (ord+1 x 1)
            F_new(1:ord+1)=1; % ������������ ��� ������ ��������� ��� ��� Newton-Raphson (ord+1 x 1)
            sigma=0.0001;
            
            while any(F_new>error); % �������� ��� loop ��� ��� ����������� ��� �� �������� ��� ������ ��������� ����� ���������� ��� �� ���� ���������
            %while sqrt(ctranspose(F_new)*F_new)>error;
                F_old=[(S-x_old(ord+1)*O)*x_old(1:ord);sum((x_old(1:ord).^2))-1]; % ������� ��������� ��� ��� Newton-Raphson - (ord+1 x 1)
                
                J_temp=[(S-x_old(ord+1)*O) -x_old(1:ord)]; % ���������� ����������� ��� ���������� ������ ��� ��� Newton-Rapshon - (ord x ord+1)
                J_old=[J_temp; 2*transpose(x_old(1:ord)) 0]; % ���������� ������� ��� ��� Newton-Rapshon - (ord+1 x ord+1)               
                
                H_old=sqrt(ctranspose(J_old)*J_old)+sqrt(sigma*diag(diag(ctranspose(J_old)*J_old)));
                %H_old=(transpose(J_old)*J_old)+(sigma*eye(ord+1));
                
                x_new=x_old-(H_old)\(transpose(J_old)*F_old); % �������� Newton-Raphson ��� ��� ������ ��� ���� x - (ord+1 x 1)
                                
                F_new=[(S-x_new(ord+1)*O)*x_new(1:ord);sum((x_new(1:ord).^2))-1]; % ����������� ��� ���� ������ ��������� ������� �� �� ��� �������� x ��� ������ ��� loop - (���� �� ������� ������ ���� ������ �� ��������������� �� ����!!)
                
                %if any(F_new>F_old)
                if sqrt(ctranspose(F_new)*F_new)>=sqrt(ctranspose(F_old)*F_old)
                    %sigma=min(10*sigma,1e16);    
                    sigma=10*sigma;
                else
                    x_old=x_new;
                    %sigma=max(0.1*sigma,1e-15);
                    sigma=0.1*sigma;
                end
            end   
            
            Tidis(:,o)=x_old(1:ord); % ���������� ��� ������������� ��������������� - (ord x ord)
            ldis(o)=x_old(ord+1); % ���������� ��� ������������� ��������� - (1 x ord)
        end
        
        %ldis=ldis.*norm; % ���������� ��� ��������������� - (1 x ord)
        ldis=(ldis+1).*norm;
        g(k,:)=sqrt(ldis); % ����������� ��� ���������� ��� propagation constant - (1 x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % ���������� ��� ������ ��������������� - (num_files x ord^2)
        end
        
    end
end

        
    
