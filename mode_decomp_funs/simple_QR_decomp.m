function [Ti,g]=simple_QR_decomp(ord,num_files,Z,Y)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα Z' για τη συγκεκριμένη συχνότητα - (ord x ord)
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα Υ' για τη συγκεκριμένη συχνότητα - (ord x ord)
    end
    
    [Tidis,ldis]=eig(Ydis*Zdis); % Υπολογισμός ιδιοδιανυσμάτων και ιδιοτιμών με βάση την QR decomposition - και οι δυο πινακες (ord x ord)
    gdis=sqrt(diag(ldis)); % Υπολογισμός του propagation constant - (1 x ord)
    
    g(k,:)=gdis; % Αποθήκευση του propagation constant - (num_files x ord)
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
    end
    
end