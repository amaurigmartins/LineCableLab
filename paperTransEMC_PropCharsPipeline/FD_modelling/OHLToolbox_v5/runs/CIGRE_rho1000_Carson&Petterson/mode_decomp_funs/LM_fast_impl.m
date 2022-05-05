function [Ti,g]=LM_fast_impl(ord,num_files,Z,Y,freq)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

O=eye(ord); % Unit matrix - (ord x ord)
error=1e-8; % Όριο σφάλματος όπου τερματίζει ο αλγόριθμος

e0=8.85418782*1e-12; % Permittivity κενού-αέρα
m0=4*pi*1e-7; % Permeability κενού-αέρα

for k=1:1:num_files % Με αύξηση της συχνότητας!
    k
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % Σύνθεση του Z' για τη παρούσα συχνότητα - (ord x ord)
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % Σύνθεση του Y' για τη παρούσα συχνότητα - (ord x ord)
    end   
    
    if k==1 % Στην πρώτη συχνότητα εφαρμόζεται QR decomposition! - Προσοχή!!! Χειροκίνητα μπορεί να εφαρμοστεί και στη δεύτερη συχνότητα κ.ο.κ.
        S=Ydis*Zdis; % Γινόμενο Y'Z' - (ord x ord)
        [Tidis,ldis]=eig(S); % QR decomposition - Όλοι οι πίνακες είναι (ord x ord)
        gdis=sqrt(diag(ldis)); % Εύρεση του propagation constant - (1 x ord)

        g(k,:)=gdis; % Αποθήκευση του propagation constant - (num_files x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
        end

    else % Στις υπόλοιπες συχνότητες εφαρμόζεται ο αλγόριθμος Newton-Raphson
        S=Ydis*Zdis; % Γινόμενο Y'Z' - (ord x ord)
        norm=-((2*pi*freq(k)).^2)*e0*m0; % Νόρμα η οποία κανονικοποιεί το γινόμενο Y'Z' ώστε να υπολογίζονται ευκολότερα οι ιδιοτιμές - Προσοχή!!! Οι ιδιοτιμές που θα υπολογιστούν χρειάζονται αντίστροφη διόρθωση με την ίδια νόρμα ενώ τα ιδιοδιανύσματα παραμένουν αναλλοίωτα αφού είναι πολλαπλασιασμός με ένα σταθερό αριθμό
        %S=S./norm; % Κανονικοποίηση του γινομένου Y'Z' - (ord x ord) 
        S=(S./norm)-O;
        
        gdis=g(k-1,:); % Εισαγωγή του propagation constant της προηγούμενης συχνότητας - Seeding!! - (1 x ord)
        ldis=gdis.^2; % Υπολογισμός των ιδιοτιμών - (1 x ord)
        %ldis=ldis./norm; % Κανονικοποίηση των ιδιοτιμών
        ldis=(ldis./norm)-1;

        for o=1:ord
            Tidis(o,:)=Ti(k-1,(o-1)*ord+1:o*ord); % Εισαγωγή του πίνακα ιδιοδιανυσμάτων της προηγούμενης συχνότητας - Seeding!! - (ord x ord)
        end
        
        for o=1:ord
            x_old=[Tidis(:,o);ldis(o)]; % Το διάνυσμα x περιέχει ένα ιδιοδιάνυσμα και την αντίστοιχη ιδιοτιμή - (ord+1 x 1)
            F_new(1:ord+1)=1; % Αρχικοποίηση του πίνακα εξισώσεων για την Newton-Raphson (ord+1 x 1)
            sigma=0.0001;
            
            while any(F_new>error); % Εκτέλεση του loop για όσο οποιοδήποτε από τα στοιχεία του πίνακα εξισώσεων είναι μεγαλύτερο από το όριο σφάλματος
            %while sqrt(ctranspose(F_new)*F_new)>error;
                F_old=[(S-x_old(ord+1)*O)*x_old(1:ord);sum((x_old(1:ord).^2))-1]; % Πίνακας εξισώσεων για την Newton-Raphson - (ord+1 x 1)
                
                J_temp=[(S-x_old(ord+1)*O) -x_old(1:ord)]; % Προσωρινός σχηματισμός του Ιακωβιανού πίνακα για την Newton-Rapshon - (ord x ord+1)
                J_old=[J_temp; 2*transpose(x_old(1:ord)) 0]; % Ιακωβιανός πίνακας για την Newton-Rapshon - (ord+1 x ord+1)               
                
                H_old=sqrt(ctranspose(J_old)*J_old)+sqrt(sigma*diag(diag(ctranspose(J_old)*J_old)));
                %H_old=(transpose(J_old)*J_old)+(sigma*eye(ord+1));
                
                x_new=x_old-(H_old)\(transpose(J_old)*F_old); % Εφαρμογή Newton-Raphson για την εύρεση του νέου x - (ord+1 x 1)
                                
                F_new=[(S-x_new(ord+1)*O)*x_new(1:ord);sum((x_new(1:ord).^2))-1]; % Υπολογισμός του νέου πίνακα εξισώσεων σύμφωνα με το νέο διάνυσμα x για έλεγχο στο loop - (Αυτό το κομμάτι κώδικα ίσως μπορεί να βελτιστοποιηθεί κι άλλο!!)
                
                %if any(F_new>F_old)
                if sqrt(ctranspose(F_new)*F_new)>=sqrt(ctranspose(F_old)*F_old)
                    %sigma=min(10*sigma,1e16);    
                    sigma=10*sigma;
                else
                    x_old=x_new;
                    %sigma=max(0.1*sigma,1e-15);
                    sigma=0.1*sigma;
                end
            end   
            
            Tidis(:,o)=x_old(1:ord); % Αποθήκευση του υπολογισμένου ιδιοδιανύσματος - (ord x ord)
            ldis(o)=x_old(ord+1); % Αποθήκευση της υπολογισμένης ιδιοτιμής - (1 x ord)
        end
        
        %ldis=ldis.*norm; % Αντιστροφή της κανονικοποίησης - (1 x ord)
        ldis=(ldis+1).*norm;
        g(k,:)=sqrt(ldis); % Υπολογισμός και αποθήκευση του propagation constant - (1 x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
        end
        
    end
end

        
    
