function [Z, Y, f, lbl] = get_ZY_from_mat(fname)

    % Load the MAT-file
    MATload = load(fname);

    % Check for the presence of required variables
    requiredVars = {'Z', 'Y', 'f'};
    for i = 1:length(requiredVars)
        if ~isfield(MATload, requiredVars{i})
            error('Missing variable: %s', requiredVars{i});
        end
    end

    if isfield(MATload,'label')
        lbl=MATload.label;
    else
        lbl='Imported from MAT-file';
    end


    % Extract variables
    Z = MATload.Z;
    Y = MATload.Y;
    f = MATload.f;

    % Validate dimensions of Z and Y
    if ndims(Z) ~= 3
        error('Variable Z must be a 3D array.');
    end
    if ndims(Y) ~= 3
        error('Variable Y must be a 3D array.');
    end
    if size(Z, 1) ~= size(Z, 2)
        error('Variable Z must have the same number of rows and columns.');
    end
    if size(Y, 1) ~= size(Y, 2)
        error('Variable Y must have the same number of rows and columns.');
    end

    % Validate dimensions consistency between Z, Y, and f
    Nph = size(Z, 1);
    NN = size(Z, 3);
    if size(Y, 1) ~= Nph || size(Y, 2) ~= Nph || size(Y, 3) ~= NN
        error('Variable Y must have the same dimensions as Z.');
    end
    if length(f) ~= NN
        error('Variable f must have the same length as the number of pages in Z and Y.');
    end

    % If all checks passed
    fprintf('Validation successful. The file %s contains valid variables.\n', fname);
end