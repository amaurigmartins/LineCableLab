function [H_mod]=calc_prop_fun(Ti,g,length,ord,num_files)

Ti_dis=zeros(ord,ord);
H_mod=zeros(num_files,ord);

for k=1:1:num_files
    
    for o=1:ord
        Ti_dis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
    end
    
    H_mod(k,:)=exp(-g(k,:)*length);
       
end
