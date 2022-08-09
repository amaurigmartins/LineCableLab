function [pol_co]=calc_poles_zeros(ord,fit)

for k=1:ord
    for o=1:ord
        [pol_co(k,o).b,pol_co(k,o).a]=residue(fit(k,o).C,fit(k,o).A,fit(k,o).D);
        pol_co(k,o).poles=roots(pol_co(k,o).a);
        pol_co(k,o).zeros=roots(pol_co(k,o).b);
    end
end
