function [V]=calc_Vnode_6ph(ord,f,Ybranch,Vo1,Vo2,Vo3,Vo4,Vo5,Vo6)

V=zeros(2*ord,max(size(f))); % (2*ord x max(size(f)))
Ybranch_dis_f=zeros(2*ord,2*ord); % (2*ord x 2*ord)
for k=1:1:max(size(f))
    
    for o=1:1:2*ord
        Ybranch_dis_f(o,:)=Ybranch(k,(o-1)*2*ord+1:o*2*ord); % Μετατροπή του πίνακα Ybranch της συγκεκριμένης συχνότητας του φάσματος σε τετραγωνικό (2*ord x 2*ord) ώστε να χρησιμοποιηθεί σε πράξεις πινάκων
    end
    
    Ired=-Ybranch_dis_f(7:length(Ybranch_dis_f),1)*Vo1(k)-Ybranch_dis_f(7:length(Ybranch_dis_f),2)*Vo2(k)-Ybranch_dis_f(7:length(Ybranch_dis_f),3)*Vo3(k)-Ybranch_dis_f(7:length(Ybranch_dis_f),4)*Vo4(k)-Ybranch_dis_f(7:length(Ybranch_dis_f),5)*Vo5(k)-Ybranch_dis_f(7:length(Ybranch_dis_f),6)*Vo6(k);
    Yred=Ybranch_dis_f(7:length(Ybranch_dis_f),7:length(Ybranch_dis_f)); % Υπολογισμός του πίνακα Υred για τη συγκεκριμένη συχνότητα φάσματος - ((2*ord)-1 x (2*ord)-1)
    %Vred=inv(Yred)*Ired; % Υπολογισμός του πίνακα Vred για τη συγκεκριμένη συχνότητα φάσματος - ((2*ord)-1 x 1) - Δηλαδή ο υπολογισμός της τάσης στους υπόλοιπους κόμβους, εκτός του γνωστού της εισόδου
    Vred=Yred\Ired;
    Vdis=[Vo1(k);Vo2(k);Vo3(k);Vo4(k);Vo5(k);Vo6(k);Vred]; % Αποθήκευση όλων των τάσεων των κόμβων για τη συγκεκριμένη συχνότητα φάσματος - (2*ord x 1)
    V(:,k)=Vdis; % Αποθήκευση όλων των τάσεων των κόμβων σε όλο το φάσμα συχνοτήτων - (2*ord x max(size(f)))
end