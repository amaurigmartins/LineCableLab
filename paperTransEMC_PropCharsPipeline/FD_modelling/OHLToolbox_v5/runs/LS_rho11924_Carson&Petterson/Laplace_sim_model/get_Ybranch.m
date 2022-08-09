function [Ybranch_dis]=get_Ybranch(ord,length,Ti_vector,Ys,Yr,gamma,Z_vector)

Z=zeros(ord,ord);
Ti=zeros(ord,ord);

for o=1:1:ord
    Z(o,:)=Z_vector(1,(o-1)*ord+1:o*ord); % Μετατροπή του πίνακα Z' της συγκεκριμένης διακριτής συχνότητας σε τετραγωνικό (ord x ord) ώστε να χρησιμοποιηθεί σε πράξεις πινάκων
    Ti(o,:)=Ti_vector(1,(o-1)*ord+1:o*ord); % Μετατροπή του πίνακα Ti της συγκεκριμένης διακριτής συχνότητας σε τετραγωνικό (ord x ord) ώστε να χρησιμοποιηθεί σε πράξεις πινάκων
end

Tv=transpose(inv(Ti)); % Υπολογισμός του αντίστοιχου πίνακα Tv για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)

%B=Tv*diag(sinh(gamma*length))*inv(diag(gamma))*inv(Tv)*Z; % Υπολογισμός του πίνακα B για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)
B=Tv*(diag(sinh(gamma*length))/(diag(gamma)))*(Tv\Z);
%B=Tv*(diag((exp(gamma*length)-exp(-gamma*length))/2)/(diag(gamma)))*(Tv\Z);

%D=Tv*diag(cosh(gamma*length))*inv(Tv); % Υπολογισμός του πίνακα D για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)
D=Tv*(diag(cosh(gamma*length))/Tv);
%D=Tv*(diag((exp(gamma*length)+exp(-gamma*length))/2)/Tv);

%Yseries=inv(B); % Υπολογισμός του πίνακα Yseries για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)
Yseries=B\(eye(ord)); %TEST
%Yshunt=inv(B)*D-inv(B); % Υπολογισμός του πίνακα Yshunt για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)
%Yshunt=(B\D)-inv(B);
Yshunt=(B\D)-Yseries; %TEST


Yself=Yseries+Yshunt; % Υπολογισμός του πίνακα Yself για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)
Ymutual=-Yseries; % Υπολογισμός του πίνακα Ymutual για τη συγκεκριμένη διακριτή συχνότητα - (ord x ord)

Yline=[Yself Ymutual;Ymutual Yself]; % Υπολογισμός του πίνακα Yline για τη συγκεκριμένη διακριτή συχνότητα - ((2*ord) x (2*ord))

Yload=[diag(Ys) zeros(ord);zeros(ord) diag(Yr)]; % Υπολογισμός του πίνακα Yload για τη συγκεκριμένη διακριτή συχνότητα - ((2*ord) x (2*ord))

Ybranch=Yline+Yload; % Υπολογισμός του πίνακα Ybranch για τη συγκεκριμένη διακριτή συχνότητα - ((2*ord) x (2*ord))

for o=1:1:2*ord
    Ybranch_dis(1,(o-1)*2*ord+1:o*2*ord)=Ybranch(o,:); % Αποθήκευση του πίνακα Ybranch για τη συγκεκριμένη διακριτή συχνότητα - (1 x (2*ord)^2) 
end