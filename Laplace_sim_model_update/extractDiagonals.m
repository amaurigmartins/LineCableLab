function diagVector = extractDiagonals(Z)
    % Get the size of the input matrix
    [N, ~, Ns] = size(Z);

    % Initialize the output matrix with Ns rows and N columns
    diagVector = zeros(Ns, N);

    % Loop over the frequency dimension (third dimension)
    for s = 1:Ns
        % Extract the diagonal elements of the s-th page
        diagElements = diag(Z(:,:,s));

        % Assign the diagonal elements to the s-th row of diagVector
        diagVector(s, :) = diagElements';
    end
end