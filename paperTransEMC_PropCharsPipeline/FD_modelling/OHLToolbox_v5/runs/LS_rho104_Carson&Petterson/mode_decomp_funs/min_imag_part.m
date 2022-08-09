function [Ti]=min_imag_part(Ti,ord,num_files)

Ti_temp=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Ti_temp(o,:)=Ti(k,(o-1)*ord+1:o*ord); % ������� ��� ������ ��������������� ��� �� ������������ ��������� - (ord x ord)
    end
    
    for o=1:ord
        Ti_real=real(Ti_temp(:,o)); % ����������� ��� ����������� ������ ��� ��������������� - (ord x 1)
        Ti_imag=imag(Ti_temp(:,o)); % ����������� ��� ����������� ������ ��� ��������������� - (ord x 1)
        
        Num=-2*sum(Ti_real.*Ti_imag); % ����������� ��� �������� ��� ��� ���������� ��� ������ ��� ������� �� ���������� ����� (scalar)
        Den=sum((Ti_real.^2)-(Ti_imag.^2)); % ����������� ��� ����������� ��� ��� ���������� ��� ������ ��� ������� �� ���������� ����� (scalar)
        
        theta=(atan(Num/Den))/2; % ����������� ��� ������ ��� ������� �� ���������� ����� ��� ��������������� (scalar)
        
        Ti_temp(:,o)=Ti_temp(:,o)*(exp(1i*theta)); % ������ ��� ��������������� ���� �� ������� �� ���������� ����� (ord x ord)
    end
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Ti_temp(o,:); % ���������� ��� ������ ��������������� - (num_files x ord^2)
    end
end