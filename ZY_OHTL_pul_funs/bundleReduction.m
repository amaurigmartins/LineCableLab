function ZR = bundleReduction(ph_order,Z)

num_ph = max(ph_order);


%% FIRST MATRIX OPERATION (Z1):

Z1 = Z;

for i = 0:num_ph
    [ph_pos] = find(ph_order == i);
    cond_per_ph = length(ph_pos);
    if ~isempty(ph_pos)
        if cond_per_ph > 1
            for j = 1:cond_per_ph - 1
                cond_col = ph_pos(1);
                subcond_col = ph_pos(j+1);

                Z1(:,subcond_col) = Z(:,subcond_col) - Z(:,cond_col);
            end
        end
    end
end

%% SECOND MATRIX OPERATION (Z1):

Z2 = Z1;

for i = 0:num_ph
    [ph_pos] = find(ph_order == i);
    cond_per_ph = length(ph_pos);
    if ~isempty(ph_pos)
        if cond_per_ph > 1
            for j = 1:cond_per_ph - 1
                cond_row = ph_pos(1);
                subcond_row = ph_pos(j+1);
                
                Z2(subcond_row,:)=Z1(subcond_row,:)-Z1(cond_row,:);
            end
        end
    end
end



nf = num_ph;
ng = size(Z,2) - num_ph;
ZC=Z2;

ZR=ZC(1:nf,1:nf)-ZC(1:nf,nf+1:nf+ng)*(ZC(nf+1:nf+ng,nf+1:nf+ng)\ZC(nf+1:nf+ng,1:nf));

end