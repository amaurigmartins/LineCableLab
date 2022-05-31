function [out] = list2sqmat(Ti,ord,freq_siz)
% Converts the lists from the toolbox into square matrices. Frequency is
% switched to represent the 1st-dimension.

for k=1:1:freq_siz
    for o=1:ord
        T_temp(o,:,k)=Ti(k,(o-1)*ord+1:o*ord);
    end
end

out=permute(T_temp,[3 1 2]);
end