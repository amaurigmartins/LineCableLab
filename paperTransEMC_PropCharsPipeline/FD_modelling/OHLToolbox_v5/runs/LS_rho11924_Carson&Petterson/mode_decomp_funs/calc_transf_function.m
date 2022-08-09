function [F_mod,F,f]=calc_transf_function(Ti,g,length,ord,freq)

f=1000:10:1000000;

g_spl=zeros(max(size(f)),ord);
Ti_spl=zeros(max(size(f)),ord^2);
Ti_dis=zeros(ord,ord);
F=zeros(max(size(f)),ord^2);

for o=1:1:ord
    a_fit(:,o)=pchip(freq,real(g(:,o)));
    b_fit(:,o)=pchip(freq,imag(g(:,o)));   
end

for o=1:1:ord^2
    Ti_real_fit(:,o)=pchip(freq,real(Ti(:,o)));
    Ti_imag_fit(:,o)=pchip(freq,imag(Ti(:,o)));
end

for o=1:1:ord
    g_spl(:,o)=ppval(a_fit(:,o),f)+1j*ppval(b_fit(:,o),f);
end

for o=1:1:ord^2  
    Ti_spl(:,o)=ppval(Ti_real_fit(:,o),f)+1j*ppval(Ti_imag_fit(:,o),f);
end

H_mod=exp(-g_spl*length);
F_mod=(2*H_mod)./(1+H_mod.^2);

for k=1:1:max(size(f))
    
    for o=1:ord
        Ti_dis(o,:)=Ti_spl(k,(o-1)*ord+1:o*ord);
    end
    
    F_dis=(transpose(inv(Ti_dis)))*(diag(F_mod(k,:)))*(transpose(Ti_dis));
    
    for o=1:ord
        F(k,(o-1)*ord+1:o*ord)=F_dis(o,:);
    end
    
end