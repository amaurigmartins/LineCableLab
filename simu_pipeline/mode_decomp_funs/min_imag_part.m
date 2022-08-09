function [Ti]=min_imag_part(Ti,ord,num_files)

Ti_temp=zeros(ord,ord); % (ord x ord)

for k=1:1:num_files
    
    for o=1:ord
        Ti_temp(o,:)=Ti(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα ιδιοδιανυσμάτων για τη συγκεκριμένη συχνότητα - (ord x ord)
    end
    
    for o=1:ord
        Ti_real=real(Ti_temp(:,o)); % Υπολογισμός του πραγματικού μέρους του ιδιοδιανύσματος - (ord x 1)
        Ti_imag=imag(Ti_temp(:,o)); % Υπολογισμός του φανταστικού μέρους του ιδιοδιανύσματος - (ord x 1)
        
        Num=-2*sum(Ti_real.*Ti_imag); % Υπολογισμός του αριθμήτη για τον υπολογισμό της γωνίας που μειώνει το φανταστικό μέρος (scalar)
        Den=sum((Ti_real.^2)-(Ti_imag.^2)); % Υπολογισμός του παρονομαστή για τον υπολογισμό της γωνίας που μειώνει το φανταστικό μέρος (scalar)
        
        theta=(atan(Num/Den))/2; % Υπολογισμός της γωνίας που μειώνει το φανταστικό μέρος του ιδιοδιανύσματος (scalar)
        
        Ti_temp(:,o)=Ti_temp(:,o)*(exp(1i*theta)); % Στροφή του ιδιοδιανύσματος ώστε να μειωθεί το φανταστικό μέρος (ord x ord)
    end
    
    for o=1:ord
        Ti(k,(o-1)*ord+1:o*ord)=Ti_temp(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
    end
end