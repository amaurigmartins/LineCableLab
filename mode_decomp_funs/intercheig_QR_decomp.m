function [Ti,g]=eigshuf_QR_decomp(ord,num_files,Z,Y)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)
YZ=zeros(ord,ord,num_files);

oldT0=[];
for k=1:1:num_files
    
    for o=1:ord
        Zdis(o,:,k)=Z(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα Z' για τη συγκεκριμένη συχνότητα - (ord x ord)
        Ydis(o,:,k)=Y(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα Υ' για τη συγκεκριμένη συχνότητα - (ord x ord)
    end
    YZ(:,:,k) = Ydis(1:ord,1:ord,k)*Zdis(1:ord,1:ord,k);
    [T0,D]=eig(YZ(:,:,k));
    [T0,D]=intercheig(T0,oldT0,D,ord,k); oldT0=T0;
    Tidis(:,:,k)=T0;
    ldis(:,:,k)=D;
end

for k=1:1:num_files
    gdis=sqrt(diag(ldis(:,:,k)));
    g(k,:)=gdis;
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:,k); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
    end
end