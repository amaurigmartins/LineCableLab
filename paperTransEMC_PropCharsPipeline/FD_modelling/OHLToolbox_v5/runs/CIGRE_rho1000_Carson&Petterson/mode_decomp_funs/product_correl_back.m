function [Ti,g]=product_correl_back(ord,num_files,Z,Y,freq)

[Ti,g]=simple_QR_decomp(ord,num_files,Z,Y); % Απλό QR decomposition χωρίς διόρθωση του eigenvector switchover - Υπολογίζει πίνακα ιδιοδιανυσμάτων Ti (current) (num_files x ord^2) και ιδιοτιμές g (sqrt) (num_files x ord)

%[Ti,g]=atp_Tsiam(Z,Y,ord,freq,num_files); % Η συνάρτηση του Τσιαμήτρου - Αποτελείται από QR decomposition και ελαχιστοποίηση του φανταστικού μέρους με βάση την μεθοδολογία του ATP/EMTP - Υπολογίζει πίνακα ιδιοδιανυσμάτων Ti (current) ((num_files x ord^2)) και ιδιοτιμές g (sqrt) (num_files x ord)

Ti_prev=zeros(ord,ord); % O πίνακας ιδιοδιανυσμάτων της παρούσας συχνότητας - Διορθωμένος απο switchover effect! - (ord x ord)
Ti_next=zeros(ord,ord); % O πίνακας ιδιοδιανυσμάτων της επόμενης συχνότητας - Μη διορθωμένος απο switchover effect! - (ord x ord)
flag=0; % Σήμανση για την ολική αποθήκευση του πίνακα ιδιοδιανυσμάτων

for k=num_files:-1:2 % Με μείωση της συχνότητας!
        
    for o=1:ord
        Ti_prev(o,:)=Ti(k,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα ιδιοδιανυσμάτων της παρούσας συχνότητας - Διορθωμένος απο switchover effect! - (ord x ord)
        Ti_next(o,:)=Ti(k-1,(o-1)*ord+1:o*ord); % Σύνθεση του πίνακα ιδιοδιανυσμάτων της επόμενης συχνότητας - Μη διορθωμένος απο switchover effect! - (ord x ord)
    end
    
    sign=1; % Σήμανση για το loop της διόρθωσης του eigenvector switchover
    while sign==1
        sign=0; % Reset της σήμανσης - Αν δεν υπάρξει κάποια αλλαγή τότε το loop θα σταματήσει
        Ti_next_mod=ctranspose(Ti_next); % Υπολογισμός του complex conjugate transpose πίνακα ιδιοδιανυσμάτων της επόμενης συχνότητας - (ord x ord)
        P=Ti_next_mod*Ti_prev; % Υπολογισμός του σχετικιστικού γινομένου - (ord x ord)
    
        for o=1:ord % Για κάθε γραμμή του σχετικιστικού πίνακα γινομένου
            [~,i]=max(P(o,:)); % Υπολογισμός της στήλης όπου βρίσκεται το μέγιστο στοιχείο της γραμμής "ο" του σχετικιστικού πίνακα γινομένου
            if i~=o % Εφόσον ο δείκτης της στήλης είναι διάφορος από το δείκτη της γραμμής
                sign=1; % Σήμανση ότι θα ξαναμπούμε στο loop της διόρθωσης αφού υπάρχει αλλαγή στα στοιχεία
                flag=1; % Σήμανση για μετέπειτα αποθήκευση αφού υπάρχει αλλαγή στα στοιχεία
                
                temp=Ti_next(:,i);
                Ti_next(:,i)=Ti_next(:,o); % Διόρθωση του switchover στον πίνακα ιδιοδιανυσμάτων
                Ti_next(:,o)=temp;
            
                temp=g(k-1,i);
                g(k-1,i)=g(k-1,o); % Διόρθωση του switchover στον πίνακα ιδιοτιμών και ολική αποθήκευση - (num_files x ord)
                g(k-1,o)=temp;
                break % Έξοδος από το for ώστε να ξαναυπολογιστεί το σχετικιστικό γινόμενο με τον εν μέρει διορθωμένο πίνακα ιδιδοδιανυσμάτων της επόμενης συχνότητας
            end
        end
    end
    
    if flag==1
        flag=0; % Reset της σήμανσης
        k-1 % Ένδειξη στο command window ότι έχουμε switchover σε αυτή τη συχνότητα
        
        for o=1:ord
            Ti(k-1,(o-1)*ord+1:o*ord)=Ti_next(o,:); % Ολική αποθήκευση του πίνακα ιδιοδιανυσμάτων αφού προηγουμένως έχουν διορθωθεί όλες οι στήλες του - (num_files x ord^2)
        end
    end
end

