function [Ti]=atp_min_imag(Ti,ord,num_files,Y)

Ydis=zeros(ord,ord); % (ord x ord)
Ti_temp=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % ������� ��� ������ Y' ��� �� ������������ ��������� - (ord x ord)
        Ti_temp(o,:)=Ti(k,(o-1)*ord+1:o*ord); % ������� ��� ������ ��������������� (current) ��� �� ������������ ��������� - (ord x ord)
    end
    
    Tv_temp=transpose(inv(Ti_temp)); % ����������� ��� ����������� ������ ��������������� ����� (voltage) (ord x ord)
    
    Y_mode_approx=(transpose(Tv_temp))*(imag(Ydis))*(Tv_temp); % ����������� ��� modal Y' �� ����� ��� ����������� ������ ��� phase Y' - �������!! ���� ����� ���������� ��� �������� ������ �� ����� ����� ��������! - (ord x ord) 
    shift=(-angle(diag(Y_mode_approx))+(pi/2))/2; % ����������� ��� ������ ��� ������� �� ���������� ����� ��� ��������������� (scalar)

    for o=1:ord
        Tv_temp(:,o)=Tv_temp(:,o)*(exp(1i*shift(o))); % % ������ ��� ��������������� ����� (voltage) ���� �� ������� �� ���������� ����� (ord x ord)
    end
    
    Ti_temp=inv(transpose(Tv_temp)); % ����������� ��� ����������� ������ ��������������� �������� (current) (ord x ord)
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Ti_temp(o,:); % ���������� ��� ������ ��������������� - (num_files x ord^2)
    end
end