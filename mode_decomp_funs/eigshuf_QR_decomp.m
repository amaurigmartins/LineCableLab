function [Ti,g]=eigshuf_QR_decomp(ord,num_files,Z,Y)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)
YZ=zeros(ord,ord,num_files);

for k=1:1:num_files
    
    for o=1:ord
        Zdis(o,:,k)=Z(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα Z' για τη συγκεκριμένη συχνότητα - (ord x ord)
        Ydis(o,:,k)=Y(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα Υ' για τη συγκεκριμένη συχνότητα - (ord x ord)
    end
    YZ(:,:,k) = Ydis(1:ord,1:ord,k)*Zdis(1:ord,1:ord,k);
end

[Tidis,ldis]=eigenshuffle(YZ);

for k=1:1:num_files
    gdis=sqrt(ldis(:,k));
    g(k,:)=gdis;
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:,k); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
    end
end