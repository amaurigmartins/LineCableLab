function [out] = list2sqmat(Ti,ord,freq_siz, reorderM)

if nargin == 3; reorderM=true; end

% Converts the lists from the toolbox into square matrices. Frequency is
% switched to represent the 1st-dimension.

if size(Ti,1)==1
    %means that the user supplied only one frequency sample, which must be
    %repeated to match freq_siz
    Ti=repmat(Ti,freq_siz,1); %repeat so we can use the same original codes from OHTL
end


for k=1:1:freq_siz
    for o=1:ord
        T_temp(o,:,k)=Ti(k,(o-1)*ord+1:o*ord);
    end
end

if reorderM
    out=permute(T_temp,[3 1 2]);
else
    out=T_temp;
end

end