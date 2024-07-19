function Z = parse_impedance_matrix(filename)
    % Read the file content
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file');
    end

    % Initialize variables
    Z = [];
    inRelevantSection = false;
    matrixSize = 0;
    currentBlock = 1;

    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        
        % Check if we've entered the relevant section
        if contains(line, 'SERIES IMPEDANCES - (ohms/kilometer)')
            inRelevantSection = true;
            continue;
        end
        
        % Check if we've exited the relevant section
        if contains(line, 'SHUNT ADMITTANCES (microsiemens/kilometer)')
            break;
        end
        
        % Skip irrelevant lines
        if ~inRelevantSection || isempty(strtrim(line)) || ...
           contains(line, 'Conductor Number:') || contains(line, 'TRALIN package')
            continue;
        end
        
        % Check for new matrix block
        if contains(line, 'MATRIX BLOCK')
            currentBlock = str2double(regexp(line, '\d+', 'match', 'once'));
            continue;
        end
        
        % Parse the line
        numbers = sscanf(line, '%f');
        
        if ~isempty(numbers)
            row = numbers(1);
            if row > matrixSize
                matrixSize = row;
            end
            for i = 1:length(numbers)-1
                col = 4 * (currentBlock - 1) + i;
                Z(row, col) = complex(numbers(i+1), 0);
            end
        end
    end

    fclose(fid);

    % Combine real and imaginary parts
    for i = 1:size(Z, 1)
        for j = 1:2:size(Z, 2)
            if j+1 <= size(Z, 2)
                Z(i,j) = complex(real(Z(i,j)), real(Z(i,j+1)));
                Z(i,j+1) = 0;
            end
        end
    end

    % Remove zero columns
    Z = Z(:, any(Z, 1));

    % Fill in the lower triangle of the matrix
    for i = 1:matrixSize
        for j = 1:i-1
            Z(i,j) = conj(Z(j,i));
        end
    end
end