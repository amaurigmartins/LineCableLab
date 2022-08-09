function [Ti,g]=SQP_calc_norm_back(ord,num_files,Z,Y,freq)

g=zeros(num_files,ord); % (num_files x ord)
Ti=zeros(num_files,ord^2); % (num_files x ord^2)
Zdis=zeros(ord,ord); % (ord x ord)
Ydis=zeros(ord,ord); % (ord x ord)

O=eye(ord); % Unit matrix - (ord x ord)

e0=8.85418782*1e-12; % Permittivity κενού-αέρα
m0=4*pi*1e-7; % Permeability κενού-αέρα

for k=num_files:-1:1 % Με μείωση της συχνότητας!
    %k % Δείκτης για να ενημερώνει στο command window σε ποια συχνότητα βρισκόμαστε
    
    for o=1:ord
        Zdis(o,:)=Z(k,(o-1)*ord+1:o*ord); % Σύνθεση του Z' για τη παρούσα συχνότητα - (ord x ord)
        Ydis(o,:)=Y(k,(o-1)*ord+1:o*ord); % Σύνθεση του Υ' για τη παρούσα συχνότητα - (ord x ord)
    end 
        
    if k==num_files || k==1 % Στην πρώτη και τελευταία συχνότητα εφαρμόζεται QR decomposition! Στην πρώτη διότι έτσι ξεκινά ο αλγόριθμος και στην τελευταία δίοτι είναι τα 0 Hz και δεν μπορεί να κανονικοποιηθεί με νόρμα το υπό εξέταση γινόμενο - Προσοχή!!! Χειροκίνητα μπορεί να εφαρμοστεί και στη δεύτερη συχνότητα κ.ο.κ.
        S=Ydis*Zdis; % Γινόμενο Υ'Ζ' - (ord x ord)
        [Tidis,ldis]=eig(S); % QR decomposition - Όλοι οι πίνακες είναι (ord x ord)
        gdis=sqrt(diag(ldis)); % Εύρεση του propagation constant - (1 x ord)

        g(k,:)=gdis; % Αποθήκευση του propagation constant - (num_files x ord)
    
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
        end

    else % Στις υπόλοιπες συχνότητες εφαρμόζεται ο αλγόριθμος SQP
        S=Ydis*Zdis; % Γινόμενο Υ'Ζ' - (ord x ord)
        norm=-((2*pi*freq(k)).^2)*e0*m0; % Νόρμα η οποία κανονικοποιεί το γινόμενο Υ'Ζ' ώστε να υπολογίζονται ευκολότερα οι ιδιοτιμές - Προσοχή!!! Οι ιδιοτιμές που θα υπολογιστούν χρειάζονται αντίστροφη διόρθωση με την ίδια νόρμα ενώ τα ιδιοδιανύσματα παραμένουν αναλλοίωτα αφού είναι πολλαπλασιασμός με ένα σταθερό αριθμό
        %S=S./norm; % Κανονικοποίηση του γινομένου Υ'Ζ' - (ord x ord) 
        S=(S./norm)-O;
        
        gdis=g(k+1,:); % Εισαγωγή του propagation constant της προηγούμενης συχνότητας - Seeding!! - (1 x ord)
        ldis=gdis.^2; % Υπολογισμός των ιδιοτιμών - (1 x ord)
        %ldis=ldis./norm; % Κανονικοποίηση των ιδιοτιμών
        ldis=(ldis./norm)-1;
        
        for o=1:ord
            Tidis(o,:)=Ti(k+1,(o-1)*ord+1:o*ord); % Εισαγωγή του πίνακα ιδιοδιανυσμάτων της προηγούμενης συχνότητας - Seeding!! - (ord x ord)
        end
                
        for o=1:ord
            x0=[real(Tidis(:,o));imag(Tidis(:,o));real(ldis(o));imag(ldis(o))]; % Το διάνυσμα x περιέχει πρώτα το πραγματικό και μετά το φανταστικό μέρος ένος ιδιοδιανύσματος, ενώ έπειτα το πραγματικό και φανταστικό μέρος της αντίστοιχη ιδιοτιμής - (2*ord+2 x 1)
            Q=[real(S) -imag(S);imag(S) real(S)]; % Βοηθητικός πίνακας για την υλοποίηση της συνάρτησης αντικειμένου (2*ord x 2*ord)
            obj=@(x)((transpose(Q*x(1:2*ord)-[x(2*ord+1)*O -x(2*ord+2)*O;x(2*ord+2)*O x(2*ord+1)*O]*x(1:2*ord)))*(Q*x(1:2*ord)-[x(2*ord+1)*O -x(2*ord+2)*O;x(2*ord+2)*O x(2*ord+1)*O]*x(1:2*ord)))/2; % Συνάρτηση αντικειμένου - Είναι βαθμωτή (scalar)! Επιστρέφει έναν πραγματικό αριθμό
            
            %c=@(x)[(2*pi*freq(k))/imag(sqrt(((x(2*ord+1)+1i*x(2*ord+2))+1)*norm));-(2*pi*freq(k))/imag(sqrt(((x(2*ord+1)+1i*x(2*ord+2)+1)*norm)));-real((sqrt((x(2*ord+1)+1i*x(2*ord+2)+1))*norm))]-[1/sqrt(e0*m0);0;0]; % Μη γραμμικοί περιορισμοί της μορφής c<=0 - (3 x 1)
            c=@(x)[];
            ceq=@(x)[sum((x(1:ord).^2)-(x(ord+1:2*ord).^2));sum(x(1:ord).*x(ord+1:2*ord))]-[1;0]; % Μη γραμμικοί περιορισμοί της μορφής ceq=0 - (2 x 1)
            
            nonlincon=@(x)deal(c(x),ceq(x)); % Ενσωμάτωση των μη γραμμικών περιορισμών σε μία συνάρτηση με μοναδικό όρισμα το διάνυσμα x
            options=optimset('Algorithm','sqp','Display','off','FinDiffType','central','MaxFunEvals',6000); % Διάφοροι παράμετροι για την επίλυση του προβλήματος της ελαχιστοποίησης της συνάρτησης-αντικειμένου
            
            x=fmincon(obj,x0,[],[],[],[],[],[],nonlincon,options); % Ελαχιστοποίηση της συνάρτησης αντικειμένου - Επιστροφή του διανύσματος x που επιτυγχάνει αυτή την ελαχιστοποίηση - (2*ord+2 x 1)
            
            Tidis(:,o)=x(1:ord)+1i*x(ord+1:2*ord); % Εξαγωγή του πίνακα ιδιοδιανυσμάτων - (ord x ord)
            ldis(o)=x(2*ord+1)+1i*x(2*ord+2); % Εξαγωγή της αντίστοιχης ιδιοτιμής - (ord x ord)
        end
        
        %ldis=ldis.*norm; % Αντιστροφή της κανονικοποίησης - (1 x ord)
        ldis=(ldis+1).*norm;
        g(k,:)=sqrt(ldis); % Υπολογισμός και αποθήκευση του propagation constant - (1 x ord)
                
        for o=1:ord
            Ti(k,(o-1)*ord+1:o*ord)=Tidis(o,:); % Αποθήκευση του πίνακα ιδιοδιανυσμάτων - (num_files x ord^2)
        end
        
    end
end
            
        