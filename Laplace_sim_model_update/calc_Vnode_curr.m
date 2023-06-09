function [V]=calc_Vnode_curr(ord,f,Ybranch,Io1)

V=zeros(2*ord,max(size(f))); % (2*ord x max(size(f)))
Ybranch_dis_f=zeros(2*ord,2*ord); % (2*ord x 2*ord)
for k=1:1:max(size(f))
    
    for o=1:1:2*ord
        Ybranch_dis_f(o,:)=Ybranch(k,(o-1)*2*ord+1:o*2*ord); % Μετατροπή του πίνακα Ybranch της συγκεκριμένης συχνότητας του φάσματος σε τετραγωνικό (2*ord x 2*ord) ώστε να χρησιμοποιηθεί σε πράξεις πινάκων
    end
    
    Ired=[Io1(k); zeros(length(Ybranch_dis_f)-1,1)];
    Vred=Ybranch_dis_f\Ired;

    V(:,k)=Vred; % Αποθήκευση όλων των τάσεων των κόμβων σε όλο το φάσμα συχνοτήτων - (2*ord x max(size(f)))
end