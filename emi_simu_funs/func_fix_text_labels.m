function [] = func_fix_text_labels(inputFilePath, outputFilePath)
% Read the input file
fid = fopen(inputFilePath, 'r');
content = fread(fid, '*char')';
fclose(fid);

% Regular expression pattern to match <text_content> tags with their inner <Text> tags
pattern = '<text_content\s+[^>]*>\s*<Text>(.*?)</Text>\s*</text_content>';

% Find matches using regular expression
matches = regexp(content, pattern, 'tokens');

% Loop through matches and replace <text_content> tags with modified content
for i = 1:length(matches)
    match = matches{i};
    replacement = match{1};  % Content within <Text> and </Text> tags
    content = strrep(content, sprintf('<Text>%s</Text>', replacement), replacement);
end

% Write the modified content to the output file
fid = fopen(outputFilePath, 'w');
fwrite(fid, content);
fclose(fid);
end