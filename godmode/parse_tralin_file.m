function [Z,Y,P] = parse_tralin_file(filename, ord)
% The ultimate parser for TRALIN files
% The filename is the path to the text file
% It will return the full NxN impedance matrix (complex values)

try
    % Open the file for reading
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Could not open the file.');
    end

    % Read the entire file into a cell array of lines
    fileContent = textscan(fileID, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fileID);

    % Extract lines from the cell array
    fileLines = fileContent{1};

    %% Impedance
    % Find the start and end lines for the impedance data
    startIdx = find(contains(fileLines, 'SERIES IMPEDANCES - (ohms/kilometer)'), 1, 'first');
    endIdx = find(contains(fileLines, 'SHUNT ADMITTANCES (microsiemens/kilometer)'), 1, 'first');

    if isempty(startIdx) || isempty(endIdx)
        error('Could not find the impedance data section in the file.');
    end

    % Extract the relevant lines
    dataLines = fileLines(startIdx:endIdx);

    % Remove useless headers and page information
    dataLines = dataLines(~contains(dataLines, {'SHUNT ADMITTANCES','SERIES IMPEDANCES','TRALIN package', 'PAGE', 'RUN', '====', '------------------------'}));
    % Remove empty lines
    dataLines = dataLines(~cellfun('isempty', strtrim(dataLines)));

    for k=1:numel(dataLines)
        dataLines{k}=removeExtraSpaces(dataLines{k});
        dataLines{k}=strrep(dataLines{k},'Conductor Number:','nan');
        dataLines{k}=strrep(dataLines{k},' +j ','+1i*');
        dataLines{k}=strrep(dataLines{k},' -j ','-1i*');
        dataLines{k}=strtrim(dataLines{k});
    end

    % Initialize variables
    block_start_indices = find(contains(dataLines, ': MATRIX BLOCK'));

    % Parse matrix blocks
    blocks = dataLines(contains(dataLines, ': MATRIX BLOCK'));
    num_blocks=numel(blocks);

    % Rearrange blocks into a fucking matrix
    Z=zeros(ord,ord);
    for j=1:num_blocks
        bs=block_start_indices(j);
        if j==num_blocks
            be=numel(dataLines);
            this_block=dataLines(bs+1:be);
        else
            be=block_start_indices(j+1);
            this_block=dataLines(bs+1:be-1);
        end

        X=[];
        for k=1:numel(this_block)
            block_split=split(this_block{k}).';
            numColumns = 5;
            % Pad with '0' if the number of elements is less than numColumns
            if numel(block_split) < numColumns
                block_split = [block_split, repmat({'0'}, 1, numColumns - numel(block_split))];
            end
            numeric_values = squeeze(cellfun(@(x) eval(x), block_split));
            X(k,:)=numeric_values;
        end
        row_indices = X(:,1); real_row_indices=(row_indices(~isnan(row_indices)));
        col_indices = X(1,:); real_col_indices=(col_indices(~isnan(col_indices)));
        X=X(2:end,2:end);
        for l=1:length(real_row_indices)
            for m=1:length(real_col_indices)
                zrow=real_row_indices(l);zcol=real_col_indices(m);
                Z(zrow,zcol)=X(l,m);
            end
        end

    end

    % Fill the lower triangular part of the matrix
    for i = 1:size(Z, 1)
        for j = i+1:size(Z, 2)
            Z(i,j) = Z(j,i);
        end
    end

    %% Admittance
    % Find the start and end lines for the impedance data
    startIdx = find(contains(fileLines, 'SHUNT ADMITTANCES (microsiemens/kilometer)'), 1, 'first');
    endIdx = find(contains(fileLines, 'SERIES ADMITTANCES (siemens.kilometer)'), 1, 'first');

    if isempty(startIdx) || isempty(endIdx)
        error('Could not find the admittance data section in the file.');
    end

    % Extract the relevant lines
    dataLines = fileLines(startIdx:endIdx);


    % Remove useless headers and page information
    dataLines = dataLines(~contains(dataLines, {'SERIES ADMITTANCES', 'SHUNT ADMITTANCES','SERIES IMPEDANCES','TRALIN package', 'PAGE', 'RUN', '====', '------------------------'}));
    % Remove empty lines
    dataLines = dataLines(~cellfun('isempty', strtrim(dataLines)));

    for k=1:numel(dataLines)
        dataLines{k}=removeExtraSpaces(dataLines{k});
        dataLines{k}=strrep(dataLines{k},'Conductor Number:','nan');
        dataLines{k}=strrep(dataLines{k},' +j ','+1i*');
        dataLines{k}=strrep(dataLines{k},' -j ','-1i*');
        dataLines{k}=strrep(dataLines{k},' +j-','-1i*');
        dataLines{k}=strtrim(dataLines{k});
    end

    % Initialize variables
    block_start_indices = find(contains(dataLines, ': MATRIX BLOCK'));

    % Parse matrix blocks
    blocks = dataLines(contains(dataLines, ': MATRIX BLOCK'));
    num_blocks=numel(blocks);

    % Rearrange blocks into a fucking matrix
    Y=zeros(ord,ord);
    for j=1:num_blocks
        bs=block_start_indices(j);
        if j==num_blocks
            be=numel(dataLines);
            this_block=dataLines(bs+1:be);
        else
            be=block_start_indices(j+1);
            this_block=dataLines(bs+1:be-1);
        end

        X=[];
        for k=1:numel(this_block)
            block_split=split(this_block{k}).';
            % Pad with '0' if the number of elements is less than numColumns
            if numel(block_split) < numColumns
                block_split = [block_split, repmat({'0'}, 1, numColumns - numel(block_split))];
            end
            numeric_values = squeeze(cellfun(@(x) eval(x), block_split));
            X(k,:)=numeric_values;
        end
        row_indices = X(:,1); real_row_indices=(row_indices(~isnan(row_indices)));
        col_indices = X(1,:); real_col_indices=(col_indices(~isnan(col_indices)));
        X=X(2:end,2:end);
        for l=1:length(real_row_indices)
            for m=1:length(real_col_indices)
                zrow=real_row_indices(l);zcol=real_col_indices(m);
                Y(zrow,zcol)=X(l,m);
            end
        end

    end

    % Fill the lower triangular part of the matrix
    for i = 1:size(Y, 1)
        for j = i+1:size(Y, 2)
            Y(i,j) = Y(j,i);
        end
    end

    %% Now the potential coefficients, not because we need, but because we can

    % Find the start and end lines for the impedance data
    startIdx = find(contains(fileLines, 'POTENTIAL COEFFICIENTS (meghoms.kilometer)'), 1, 'first');
    endIdx = find(contains(fileLines, 'SERIES IMPEDANCES - (ohms/kilometer)'), 1, 'first');

    if isempty(startIdx) || isempty(endIdx)
        error('Could not find the admittance data section in the file.');
    end

    % Extract the relevant lines
    dataLines = fileLines(startIdx:endIdx);


    % Remove useless headers and page information
    dataLines = dataLines(~contains(dataLines, {'POTENTIAL COEFFICIENTS', 'SERIES ADMITTANCES', 'SHUNT ADMITTANCES','SERIES IMPEDANCES','TRALIN package', 'PAGE', 'RUN', '====', '------------------------'}));
    % Remove empty lines
    dataLines = dataLines(~cellfun('isempty', strtrim(dataLines)));

    for k=1:numel(dataLines)
        dataLines{k}=removeExtraSpaces(dataLines{k});
        dataLines{k}=strrep(dataLines{k},'Conductor Number:','nan');
        dataLines{k}=strrep(dataLines{k},' +j ','+1i*');
        dataLines{k}=strrep(dataLines{k},' -j ','-1i*');
        dataLines{k}=strrep(dataLines{k},' +j-','-1i*');
        dataLines{k}=strtrim(dataLines{k});
    end

    % Initialize variables
    block_start_indices = find(contains(dataLines, ': MATRIX BLOCK'));

    % Parse matrix blocks
    blocks = dataLines(contains(dataLines, ': MATRIX BLOCK'));
    num_blocks=numel(blocks);

    % Rearrange blocks into a fucking matrix
    P=zeros(ord,ord);
    for j=1:num_blocks
        bs=block_start_indices(j);
        if j==num_blocks
            be=numel(dataLines);
            this_block=dataLines(bs+1:be);
        else
            be=block_start_indices(j+1);
            this_block=dataLines(bs+1:be-1);
        end

        X=[];
        for k=1:numel(this_block)
            block_split=split(this_block{k}).';
            % Pad with '0' if the number of elements is less than numColumns
            if numel(block_split) < numColumns
                block_split = [block_split, repmat({'0'}, 1, numColumns - numel(block_split))];
            end
            numeric_values = squeeze(cellfun(@(x) eval(x), block_split));
            X(k,:)=numeric_values;
        end
        row_indices = X(:,1); real_row_indices=(row_indices(~isnan(row_indices)));
        col_indices = X(1,:); real_col_indices=(col_indices(~isnan(col_indices)));
        X=X(2:end,2:end);
        for l=1:length(real_row_indices)
            for m=1:length(real_col_indices)
                zrow=real_row_indices(l);zcol=real_col_indices(m);
                P(zrow,zcol)=X(l,m);
            end
        end

    end

    % Fill the lower triangular part of the matrix
    for i = 1:size(Y, 1)
        for j = i+1:size(P, 2)
            P(i,j) = P(j,i);
        end
    end

    %% Let us return decent units because we are not complete retards
    Z = Z/1000;
    Y = Y*1e-6/1000;
    P = P*1e6*1000; %TESTME
catch ME
    fprintf('Error occurred: %s\n', ME.message);
    Z = [];
    Y = [];
    P = [];
end

end
