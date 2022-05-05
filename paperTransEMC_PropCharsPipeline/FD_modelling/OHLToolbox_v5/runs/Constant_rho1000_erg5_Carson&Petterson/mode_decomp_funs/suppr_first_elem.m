function [Ti]=suppr_first_elem(Ti,ord,num_files)

Ti_temp=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Ti_temp(o,:)=Ti(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα ιδιοδιανυσμάτων για τη συγκεκριμένη συχνότητα - (ord x ord)
    end
    
    suppr=Ti_temp(1,:); % Εκλογή της πρώτης γραμμής του πίνακα ιδιοδιανυσμάτων ως αναφορά για τη διαίρεση - (1 x ord)
    for o=1:ord
        Ti_temp(o,:)=Ti_temp(o,:)./suppr; % Διαίρεση όλων των γραμμών του πίνακα μετασχηματισμού με την συγκεκριμένη γραμμή - (ord x ord)
    end
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Ti_temp(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
    end
end