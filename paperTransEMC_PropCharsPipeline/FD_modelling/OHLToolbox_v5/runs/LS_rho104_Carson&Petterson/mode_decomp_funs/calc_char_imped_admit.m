function [Zch_mod,Ych_mod,Zch,Ych]=calc_char_imped_admit(Ti,Z,Y,ord,num_files)

Ti_dis=zeros(ord,ord);
Z_dis=zeros(ord,ord);
Y_dis=zeros(ord,ord);
Zch_mod=zeros(num_files,ord);
Ych_mod=zeros(num_files,ord);
Zch=zeros(num_files,ord^2);
Ych=zeros(num_files,ord^2);

for k=1:1:num_files
    
    for o=1:ord
        Z_dis(o,:)=Z(k,(o-1)*ord+1:o*ord);
        Y_dis(o,:)=Y(k,(o-1)*ord+1:o*ord);
        Ti_dis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
    end
    
    Zmodal=(transpose(Ti_dis))*(Z_dis)*(Ti_dis);
    Ymodal=(inv(Ti_dis))*(Y_dis)*(transpose(inv(Ti_dis)));
    
    Zch_mod_dis=transpose(diag(sqrt(Zmodal)./sqrt(Ymodal)));
    Ych_mod_dis=transpose(diag(sqrt(Ymodal)./sqrt(Zmodal)));
        
    Zch_dis=(transpose(inv(Ti_dis)))*(diag(Zch_mod_dis))*(inv(Ti_dis));
    Ych_dis=(Ti_dis)*(diag(Ych_mod_dis))*(transpose(Ti_dis));
    
    Zch_mod(k,:)=Zch_mod_dis;
    Ych_mod(k,:)=Ych_mod_dis;
    
    for o=1:ord
        Zch(k,(o-1)*ord+1:o*ord)=Zch_dis(o,:);
        Ych(k,(o-1)*ord+1:o*ord)=Ych_dis(o,:);
    end
    
end