%V_S=[exp(1i*(0)); exp(1i*(-2*pi/3)); exp(1i*(2*pi/3));];
%V_S=[1; 1; 1;];
V_S=[1; 0; 0;];

F_dis=zeros(ord,ord);
V_R=zeros(ord,max(size(f)));

for k=1:1:max(size(f))
    
    for o=1:ord
        F_dis(o,:)=F(k,(o-1)*ord+1:o*ord);
    end
    
    V_R_dis=F_dis*V_S;
    
    V_R(:,k)=V_R_dis;
    
end

figure(1)
for k=1:ord
    subplot(2,1,1)
    for o=1:ord
        semilogx(f,abs(V_R(o,:)))
        hold all
    end
    hold off
    subplot(2,1,2)
    for o=1:ord
        semilogx(f,radtodeg(unwrap(angle(V_R(o,:)))))
        hold all
    end
    hold off
end