function [V_mod]=calc_modal_variables(ord,f,Ti,V)

Ti_dis=zeros(ord,ord);
V_mod=zeros(ord,max(size(f)));

for k=1:1:max(size(f))
    
    for o=1:1:ord
        Ti_dis(o,:)=Ti(k,(o-1)*ord+1:o*ord);
    end
    
    V_mod(1:ord,k)=transpose(Ti_dis)*V(ord+1:2*ord,k);
end
