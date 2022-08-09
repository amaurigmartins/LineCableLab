function [Ybranch_dis]=get_Ybranch_scaled(ord,length,Ti_vector,Ys,Yr,gamma,Z_vector)

Z=zeros(ord,ord);
Ti=zeros(ord,ord);

for o=1:1:ord
    Z(o,:)=Z_vector(1,(o-1)*ord+1:o*ord); % Μετατροπή του πίνακα Z' της συγκεκριμένης διακριτής συχνότητας σε τετραγωνικό (ord x ord) ώστε να χρησιμοποιηθεί σε πράξεις πινάκων
    Ti(o,:)=Ti_vector(1,(o-1)*ord+1:o*ord); % Μετατροπή του πίνακα Ti της συγκεκριμένης διακριτής συχνότητας σε τετραγωνικό (ord x ord) ώστε να χρησιμοποιηθεί σε πράξεις πινάκων
end

Tv=transpose(inv(Ti)); % Υπολογισμός του αντίστοιχου πίνακα Tv για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)

h=exp(-gamma*length);

alpha=(gamma.*(1+h.^2))./(1-h.^2);
beta=(-2*gamma.*h)./(1-h.^2);

%A=inv(Z)*Tv*diag(alpha)*inv(Tv);
A=(Z\eye(ord))*Tv*diag(alpha)*(Tv\eye(ord));

%B=inv(Z)*Tv*diag(beta)*inv(Tv);
B=(Z\eye(ord))*Tv*diag(beta)*(Tv\eye(ord));

Yline=[A B;B A]; % Υπολογισμός του πίνακα Yline για τη συγκεκριμένη διακριτή συχνότητα - ((2*ord) x (2*ord))

Yload=[diag(Ys) zeros(ord);zeros(ord) diag(Yr)]; % Υπολογισμός του πίνακα Yload για τη συγκεκριμένη διακριτή συχνότητα - ((2*ord) x (2*ord))

Ybranch=Yline+Yload; % Υπολογισμός του πίνακα Ybranch για τη συγκεκριμένη διακριτή συχνότητα - ((2*ord) x (2*ord))

for o=1:1:2*ord
    Ybranch_dis(1,(o-1)*2*ord+1:o*2*ord)=Ybranch(o,:); % Αποθήκευση του πίνακα Ybranch για τη συγκεκριμένη διακριτή συχνότητα - (1 x (2*ord)^2) 
end