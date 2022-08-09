function [Ti]=atp_min_imag(Ti,ord,num_files,Y)

Ydis=zeros(ord,ord); % (ord x ord)
Ti_temp=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα Y' για τη συγκεκριμένη συχνότητα - (ord x ord)
        Ti_temp(o,:)=Ti(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα ιδιοδιανυσμάτων (current) για τη συγκεκριμένη συχνότητα - (ord x ord)
    end
    
    Tv_temp=transpose(inv(Ti_temp)); % Υπολογισμός του αντίστοιχου πίνακα ιδιοδιανυσμάτων τάσης (voltage) (ord x ord)
    
    Y_mode_approx=(transpose(Tv_temp))*(imag(Ydis))*(Tv_temp); % Υπολογισμός του modal Y' με χρήση του φανταστικού μέρους του phase Y' - Προσοχή!! Αυτό είναι προσέγγιση και περιέχει σφάλμα το οποίο είναι αποδεκτό! - (ord x ord) 
    shift=(-angle(diag(Y_mode_approx))+(pi/2))/2; % Υπολογισμός της γωνίας που μειώνει το φανταστικό μέρος του ιδιοδιανύσματος (scalar)

    for o=1:ord
        Tv_temp(:,o)=Tv_temp(:,o)*(exp(1i*shift(o))); % % Στροφή του ιδιοδιανύσματος τάσης (voltage) ώστε να μειωθεί το φανταστικό μέρος (ord x ord)
    end
    
    Ti_temp=inv(transpose(Tv_temp)); % Υπολογισμός του αντίστοιχου πίνακα ιδιοδιανυσματων ρεύματος (current) (ord x ord)
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Ti_temp(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
    end
end