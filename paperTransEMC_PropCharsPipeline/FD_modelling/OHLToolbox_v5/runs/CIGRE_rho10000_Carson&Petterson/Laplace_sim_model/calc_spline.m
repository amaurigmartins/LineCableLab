function [Ybranch]=calc_spline(ord,freq,f,Ybranch_dis_tot)
    
for o=1:1:(2*ord)^2
   Ybranch_real_fit(:,o)=spline(freq,real(Ybranch_dis_tot(:,o))); % Spline στο πραγματικό μέρος κάθε στήλης του πίνακα Ybranch_dis_tot - Δεν ορίζεται διάσταση, είναι τύπου struct
   Ybranch_imag_fit(:,o)=spline(freq,imag(Ybranch_dis_tot(:,o))); % Spline στο φανταστικό μέρος κάθε στήλης του πίνακα Ybranch_dis_tot - Δεν ορίζεται διάσταση, είναι τύπου struct 
end

Ybranch=zeros(max(size(f)),(2*ord)^2); % (max(size(f)) x (2*ord)^2)
for o=1:1:(2*ord)^2
    Ybranch(:,o)=ppval(Ybranch_real_fit(:,o),f)+1j*ppval(Ybranch_imag_fit(:,o),f); % Ppval στις Spline για να γίνει η αριθμητική προσέγγιση σε όλο το φάσμα συχνοτήτων - (max(size(f)) x (2*ord)^2)
end