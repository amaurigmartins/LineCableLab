Zdis=zeros(ord,ord);
Ydis=zeros(ord,ord);
S=zeros(num_files,ord^2);

e0=8.85418782*1e-12; % Permittivity κενού-αέρα
m0=4*pi*1e-7; % Permeability κενού-αέρα
O=eye(ord);

for k=1:1:num_files
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord);
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord);
    end
    
    Sdis=Ydis*Zdis;
    norm=-((2*pi*freq(k)).^2)*e0*m0;
    
    Sdis=(Sdis./norm)-O;
    
    for o=1:ord
        S(k,(o-1)*ord+1:o*ord)=Sdis(o,:); 
    end
end

figure(1)
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,real(S(:,k+(o-1)*ord)))
        hold all
    end
    hold off
    subplot(ord,2,2*k)
    for o=1:ord
        semilogx(freq,imag(S(:,k+(o-1)*ord)))
        hold all
    end
    hold off
end

figure(2)
for k=1:ord
    subplot(ord,2,2*k-1)
    for o=1:ord
        semilogx(freq,abs(S(:,k+(o-1)*ord)))
        hold all
    end
    hold off
    subplot(ord,2,2*k)
    for o=1:ord
        semilogx(freq,radtodeg(unwrap(angle(S(:,k+(o-1)*ord)))))
        hold all
    end
    hold off
end