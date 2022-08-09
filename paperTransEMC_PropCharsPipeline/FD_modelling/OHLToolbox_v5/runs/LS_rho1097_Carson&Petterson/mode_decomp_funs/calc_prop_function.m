function [H_mod,P_mod,H,P,delay]=calc_prop_function(Ti,g,length,ord,num_files,freq)

Ti_dis=zeros(ord,ord);
H_mod=zeros(num_files,ord);
P_mod=zeros(num_files,ord);
P=zeros(num_files,ord^2);
H=zeros(num_files,ord^2);

vel(1:ord)=(2*pi*freq(num_files))./imag(g(num_files,:));
delay(1:ord)=length./vel;

for k=1:1:num_files
    
    for o=1:ord
        Ti_dis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
    end
    
    H_mod(k,:)=exp(-g(k,:)*length);
    %H_mod(k,:)=abs(exp(-g(k,:)*length)).*exp(1i*(-imag(g(k,:)*length)));
    P_mod(k,:)=H_mod(k,:).*exp(1i*2*pi*freq(k).*delay);
    
    P_dis=(Ti_dis)*(diag(P_mod(k,:)))*(inv(Ti_dis));
    %H_dis=(Ti_dis)*(diag(H_mod(k,:)))*(inv(Ti_dis));
    H_dis=(Ti_dis)*(diag(exp(-1i*2*pi*freq(k).*delay)))*(inv(Ti_dis))*P_dis;
    
    for o=1:ord
        P(k,(o-1)*ord+1:o*ord)=P_dis(o,:);
        H(k,(o-1)*ord+1:o*ord)=H_dis(o,:);
    end
        
end