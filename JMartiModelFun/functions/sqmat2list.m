function [out] = sqmat2list(Ti,ord,freq_siz)
% Converts square matrices into lists usable by the toolbox. Each line
% represents a frequency sample.

if ndims(Ti)==2
    %means that the user supplied only one frequency sample, which must be
    %repeated to match freq_siz
    Ti=repmat(Ti,1,1,freq_siz); %repeat so we can use the same original codes from OHTL
    Ti=permute(Ti,[3 1 2]); %make frequency the first dimension
end

for k=1:1:freq_siz
    for o=1:ord
        T_temp(k,(o-1)*ord+1:o*ord)=Ti(k,o,:); % ÁðïèÞêåõóç ôïõ ðßíáêá éäéïäéáíõóìÜôùí - (num_files x ord^2)
    end
end

out=T_temp;
end