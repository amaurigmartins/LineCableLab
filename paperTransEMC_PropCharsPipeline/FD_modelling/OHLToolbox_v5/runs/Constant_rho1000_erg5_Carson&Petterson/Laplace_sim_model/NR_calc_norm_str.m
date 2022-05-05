function [Ti,g]=NR_calc_norm_str(ord,num_files,Z,Y,freq)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

O=eye(ord); % Unit matrix - (ord x ord)
error=1e-8; % Όριο σφάλματος όπου τερματίζει ο αλγόριθμος

e0=8.85418782*1e-12; % Permittivity κενού-αέρα
m0=4*pi*1e-7; % Permeability κενού-αέρα

for k=1:1:num_files % Με αύξηση της συχνότητας!
    
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % Σύνθεση του Z' για τη παρούσα συχνότητα - (ord x ord)
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % Σύνθεση του Y' για τη παρούσα συχνότητα - (ord x ord)
    end   
    
    if k==1 || k==2 % Στην πρώτη συχνότητα εφαρμόζεται QR decomposition! - Προσοχή!!! Χειροκίνητα μπορεί να εφαρμοστεί και στη δεύτερη συχνότητα κ.ο.κ.
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
            x=[Tidis(:,o);ldis(o)]; % Το διάνυσμα x περιέχει ένα ιδιοδιάνυσμα και την αντίστοιχη ιδιοτιμή - (ord+1 x 1)
            F(1:ord+1)=1; % Αρχικοποίηση του πίνακα εξισώσεων για την Newton-Raphson (ord+1 x 1)
            
            while any(F>error); % Εκτέλεση του loop για όσο οποιοδήποτε από τα στοιχεία του πίνακα εξισώσεων είναι μεγαλύτερο από το όριο σφάλματος
                F=[(S-x(ord+1)*O)*x(1:ord);sum((x(1:ord).^2))-1]; % Πίνακας εξισώσεων για την Newton-Raphson - (ord+1 x 1)
                J_temp=[(S-x(ord+1)*O) -x(1:ord)]; % Προσωρινός σχηματισμός του Ιακωβιανού πίνακα για την Newton-Rapshon - (ord x ord+1)
                J=[J_temp; 2*transpose(x(1:ord)) 0]; % Ιακωβιανός πίνακας για την Newton-Rapshon - (ord+1 x ord+1)
                x=x-(J\F); % Εφαρμογή Newton-Raphson για την εύρεση του νέου x - (ord+1 x 1)
                F=[(S-x(ord+1)*O)*x(1:ord);sum((x(1:ord).^2))-1]; % Υπολογισμός του νέου πίνακα εξισώσεων σύμφωνα με το νέο διάνυσμα x για έλεγχο στο loop - (Αυτό το κομμάτι κώδικα ίσως μπορεί να βελτιστοποιηθεί κι άλλο!!)
            end   
            
            Tidis(:,o)=x(1:ord); % Αποθήκευση του υπολογισμένου ιδιοδιανύσματος - (ord x ord)
            ldis(o)=x(ord+1); % Αποθήκευση της υπολογισμένης ιδιοτιμής - (1 x ord)
        end
        
        %ldis=ldis.*norm; % Αντιστροφή της κανονικοποίησης - (1 x ord)
        ldis=(ldis+1).*norm;
        g(k,:)=sqrt(ldis); % Υπολογισμός και αποθήκευση του propagation constant - (1 x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
        end
        
    end
end


        
    
