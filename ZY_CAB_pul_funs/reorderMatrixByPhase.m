function [Z_reordered, ph_reordered] = reorderMatrixByPhase(Z, ph_order)
    % Initialize an empty array to store the indices of the reordered rows/columns
    reordered_indices = [];
    
    % Get the unique phases from ph_order
    unique_phases = unique(ph_order(ph_order > 0)); % Ignore phase 0 for now
    
    % First, we process one row for each unique phase
    for phase = unique_phases'
        % Get the first row index corresponding to the current phase
        idx = find(ph_order == phase, 1, 'first');
        if ~isempty(idx)
            reordered_indices = [reordered_indices, idx];
        end
    end
    
    % Now, append all remaining rows of each phase
    for phase = unique_phases'
        % Get all row indices corresponding to the current phase, excluding the first one
        idx = find(ph_order == phase);
        idx(1) = [];  % Remove the first row (it was already added above)
        reordered_indices = [reordered_indices, idx'];
    end
    
    % Finally, append all rows/columns corresponding to phase 0 at the end
    zero_phase_indices = find(ph_order == 0);
    reordered_indices = [reordered_indices, zero_phase_indices'];
    
    % Now reorder both rows and columns of Z based on the reordered indices
    Z_reordered = Z(reordered_indices, reordered_indices);
    
    % Also reorder the phase order based on the same indices
    ph_reordered = ph_order(reordered_indices);
end
