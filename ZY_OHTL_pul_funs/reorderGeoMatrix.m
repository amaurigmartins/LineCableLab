function MF = reorderGeoMatrix(M)

M = sortrows(M,1);

zero_ph_row = find(M(:,1) == 0);
if ~isempty(zero_ph_row)
     M0 = M(zero_ph_row(1):zero_ph_row(end),:);
     M1 = M(zero_ph_row(end)+1:end,:);
     PH = unique(M(zero_ph_row(end)+1:end,1));
     
else
    M1 = M;
    M0 = [];
end

MF = [];
k = 1;
i = 1;
N = size(M1,1);

Nph = unique(M1(:,1));
Nph = size(Nph(Nph~=0),1);

for i = 1:2
    for j = 1:Nph
        if i == 1
            ph_row = find(M1(:,1)==j);
            if ~isempty(ph_row)
                MF(j,:) = M1(ph_row(1),:);
                M1(ph_row(1),:) = [];
            end
        else
            ph_row = find(M1(:,1)==j);
            if ~isempty(ph_row)
                MF_aux = M1(ph_row,:);
                M1(ph_row,:) = [];
                MF = [MF;MF_aux];
            end
        end
    end
end
        
    

% while size(MF,1) ~= N
%    flag = 0;
%     while flag == 0
%        ph_row = find(M1(:,1)==k);
%         if ~isempty(ph_row)
%             MF(i,:) = M1(ph_row(1),:);
%             flag = 1;
%             k = k + 1;
%             i = i + 1; 
%             M1(ph_row(1),:) = [];
%         else
%             k = min(M1(:,1));
%             flag = 1;
%         end
%     end 
% end
    
MF = [MF; M0];

end
