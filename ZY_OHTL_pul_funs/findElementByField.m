function index = findElementByField(structArray, fieldName, testString)
%findStructElementByField - Finds the index of the struct element whose field matches the test string
%
% Syntax: index = findStructElementByField(structArray, fieldName, testString)
%
% Inputs:
%    structArray - An array of structs
%    fieldName - The name of the field to test
%    testString - The string to match against the field's value
%
% Outputs:
%    index - The index of the matching struct element. Returns the next value in array if no match is found.

    % Check input validity
    if ~isstruct(structArray)
        error('Input must be an array of structs.');
    end

    % Use arrayfun to check each element
    index = find(arrayfun(@(x) strcmp(x.(fieldName), testString), structArray), 1);

    if isempty(index); index=numel(structArray)+1; end
end
