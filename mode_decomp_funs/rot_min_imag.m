function [Ti]=rot_min_imag(Ti,ord,num_files)

Ti_temp=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Ti_temp(o,:)=Ti(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα ιδιοδιανυσμάτων (current) για τη συγκεκριμένη συχνότητα - (ord x ord)
    end
    
    Ti_temp=rot(Ti_temp);
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Ti_temp(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
    end
end