% ParametricRunBuilder.m
%
% DESCRIPTION (yeah, we still document things like civilized barbarians):
%   Implements a  kludgy-but-effective parametric runner with:
%   - Ranged parameters (cell arrays of strings) producing full Cartesian combos
%   - Deterministic hash-like job IDs: BatchID + '_' + raw values ('.' -> '_')
%   - Non-destructive workflow: copies original job scripts to RunBatchJob_<hashID>.m
%   - Uses function override_params() to inject per-run singletons + jobid
%   - Generates a RunAllJobs.m that executes every generated job
%
% HOW TO USE:
%   1) Configure "CONFIGURATION" below:
%       - BatchID (prefix)
%       - var_ranges (variables to change in the original files, values can be single strings or cellstr ranges)
%       - (optional) static_vars to slam constants across all runs
%       - srcdir, tgtdir, fname (or let dialogs annoy you)
%   2) Run this script. It spawns one job file per combo per original, and runs them.
%
% Pro-tip:
%   Values in var_ranges must be strings of MATLAB code.
%   For string literals you’re writing into files, remember the joy of doubled quotes.
%   Example: {'jobtitle', '''oh_look_a_string'''}  % two single-quotes inside to emit one.
%
% Workflow:
%
% var_ranges  ──┐
%               │  expand_var_ranges  →  [combo1, combo2, ... comboN]
% static_vars ──┘                                  │
%                          for each original job   │  for each combo
%                                   ┌──────────────┴──────────────┐
%                                   │ copy → RunBatchJob_<hash>.m │
%                                   │ override_params on the copy │
%                                   │ append run stanza to master │
%                                   └──────────────┬──────────────┘
%                                                  ↓
%                                          RunAllJobs.m
%                                                  ↓
%                              runs → copies *_NLTSimu.mat, *_LineParam.mat, plots


%% OPTIONAL - MAP FOLDER TO VIRTUAL DRIVE ------------------------------------------------------------
driveLetter = 'Z';
folderPath = ''; % keep empty unless you're into Windows cosplay

if ~isempty(folderPath) && isfolder(folderPath)
    fprintf('Attempting to map drive %s: to path: %s\n', driveLetter, folderPath);
    [~, ~] = mapDrive(folderPath, driveLetter);
else
    fprintf('Drive mapping skipped: folderPath is empty or not a valid folder.\n');
end

%% 1. INITIAL WORKSPACE CLEANUP ----------------------------------------------------------------------
close all;
clearvars -except folderPath driveLetter; % stop obliterating everything like a toddler
delete(allchild(0));
clc;

%% 2. CONFIGURATION ----------------------------------------------------------------------------------
BatchID = 'testBatch'; % prefix that will be stitched to each combo value to mint hashID

% --- Ranged & singleton parameters ------------------------------------------------------------------
% Values are STRINGS OF MATLAB CODE that will be written as the right-hand side.
% If you want a literal string in the target file, you need doubled single-quotes.
% Example: '''abc''' emits 'abc' in the file.
var_ranges = {
    'runLineParametersCalc', 'logical(1)';
    'runSimuTool',           'logical(1)';
    'ZYprnt',                'logical(1)';
    'Modprnt',               'logical(1)';
    'f',                     'transpose(logspaceppd(-1,6,20))';
    'ZYsave',                'logical(1)';
    'Figsave',               'logical(1)';
    'export2EMTP',           'logical(0)';
    'export2XML',            'logical(0)';
    'export2ULM',            'logical(0)';
    'decmp_flag',            {'9','10'};              % <---- example range
    'rho',                   {'100','200'};           % <---- another range
};

% If you want some constants slammed on every run regardless of ranges, do it here.
% (Still strings of MATLAB code, same convention.)
static_vars = {
    % 'someOtherFlag', '42';
};

% --- Directory & file pattern -----------------------------------------------------------------------
% srcdir = 'C:\path\to\source';     % leave undefined to get UI
% tgtdir = 'C:\path\to\target';     % leave undefined to get UI
% fname  = 'RunJob_*1000*.m';       % leave undefined to get UI
post_process = '';                  % leave empty if not used

% nukes target output .mat and plots before starting — your circus, your fire
force_cleanup = true;

%% 3. SETUP & VALIDATION -----------------------------------------------------------------------------
if ~exist('srcdir','var')
    srcdir = uigetdir('', 'Select source directory for recursive search');
    if isequal(srcdir,0); error('No directory selected. Exiting...'); end
end
cd(srcdir);

if ~exist('tgtdir','var')
    tgtdir = uigetdir('', 'Select target directory for copied files');
    if isequal(tgtdir,0); error('No target directory selected. Exiting...'); end
end

if ~exist('post_process','var')
    post_process = getFilePath(); % your existing helper
end

if ~exist('fname','var')
    fname = getStringInput();
    if isempty(fname); fname = '*.m'; end
end

if force_cleanup
    cleanup_directory(tgtdir); % your existing helper — deletes .mat + plots
end

%% 4. FIND JOB FILES ---------------------------------------------------------------------------------
fprintf('Searching for files matching "%s" in "%s"...\n', fname, srcdir);
files = dir(fullfile(srcdir, '**', fname));
if isempty(files)
    warning('No matching files found in %s. Exiting.', srcdir);
    return;
end
fprintf('Found %d job files. Proceeding to spawn parametric abominations...\n', numel(files));

%% 5. EXPAND RANGES TO COMBINATIONS ------------------------------------------------------------------
% This gives:
%   combKeys  : Nx1 cell of cell arrays (each inner cell is {varName, chosenValueString})
%   suffixes  : Nx1 cell of strings used for hash tail based on raw chosen values ('.'->'_')
%   flatNames : 1xK cell of var names (order preserved)
%   flatSets  : 1xK cell of 1xMi cell arrays of candidate value strings
[combKeys, suffixes, flatNames, flatSets] = expand_var_ranges(var_ranges);

fprintf('Total combinations from var_ranges: %d\n', numel(combKeys));

%% 6. GENERATE THE MASTER EXECUTION SCRIPT -----------------------------------------------------------
scriptName = 'RunAllJobs.m';
fid = fopen(scriptName, 'w');
if fid == -1
    error('Cannot create master execution script. Grow opposable thumbs and check permissions.');
end

fprintf(fid, '%% Auto-generated script to run all jobs sequentially\n');
fprintf(fid, '%% Generated by ParametricRunBuilder.m on %s\n\n', datestr(now));
fprintf(fid, 'warning(''on'');\n');

generatedCount = 0;

%% 7. FOR EACH ORIGINAL JOB, FOR EACH COMBO: COPY, OVERRIDE, QUEUE -----------------------------------
for ff = 1:numel(files)
    origPath = fullfile(files(ff).folder, files(ff).name);
    for cc = 1:numel(combKeys)
        % Build hashID = BatchID + concatenated '_' + raw values ('.'->'_')
        hashTail = suffixes{cc};
        hashID   = [BatchID, hashTail];

        % forge target job file name next to the original
        outName  = ['RunBatchJob_', hashID, '.m'];
        outPath  = fullfile(files(ff).folder, outName);

        % 1) Copy original → new job script (non-destructive, unlike your average coder)
        copy_job_script(origPath, outPath);

        % 2) Build newvars by flattening the combo + static_vars + injected jobid
        % combo entries are already strings of code; jobid must be a STRING LITERAL in the file
        % so we wrap in doubled quotes inside the line content.
        comboPairs = combKeys{cc}; % cell Nx2: {varName, rhsString}
        newvars = [comboPairs; static_vars; {'jobid', ['''' hashID '''']}];

        % 3) Apply your existing line-smasher
        override_params(outPath, newvars);

        % 4) Append execution stanza to master script
        fprintf(fid, '%% --- Job %d/%d | Combo %d/%d: %s --> %s ---\n', ...
            ff, numel(files), cc, numel(combKeys), files(ff).name, outName);
        fprintf(fid, 'cd(''%s'');\n', files(ff).folder);
        fprintf(fid, 'try\n');
        fprintf(fid, '    fprintf(''\\n\\n======================================\\n'');\n');
        fprintf(fid, '    fprintf(''Now running: %s (combo %d/%d)\\n'');\n', outName, cc, numel(combKeys));
        fprintf(fid, '    fprintf(''======================================\\n\\n'');\n');
        fprintf(fid, '    run(''%s'');\n', outName);
        fprintf(fid, '    if exist(''jobid'',''var'')\n');
        fprintf(fid, '        copyfile([jobid ''_NLTSimu.mat''], ''%s'', ''f'');\n', tgtdir);
        fprintf(fid, '        fprintf(''Copied %%s_NLTSimu.mat to target directory.\\n'', jobid);\n');
        fprintf(fid, '        copyfile([jobid ''_LineParam.mat''], ''%s'', ''f'');\n', tgtdir);
        fprintf(fid, '        fprintf(''Copied %%s_LineParam.mat to target directory.\\n'', jobid);\n');
        fprintf(fid, '        copyfile(''./plots'', fullfile(''%s'',''plots''));\n', tgtdir);
        fprintf(fid, '        fprintf(''Copied plots subfolder to target directory.\\n'');\n');
        fprintf(fid, '    else\n');
        fprintf(fid, '        warning(''jobid variable not defined after running %s. Results may not have been copied.'');\n', outName);
        fprintf(fid, '    end\n');
        fprintf(fid, 'catch ME\n');
        fprintf(fid, '    warning(''An error occurred in job %s: %%s'', ME.message);\n', outName);
        fprintf(fid, 'end\n\n');

        generatedCount = generatedCount + 1;
    end
end

fprintf(fid, 'fprintf(''\\n\\nAll jobs completed!\\n'');\n');

if ~isempty(post_process) && isfile(post_process)
    fprintf('Adding post-processing step...\n');
    fprintf(fid, 'fprintf(''Running post-processing script...\\n'');\n');
    fprintf(fid, 'batch_tag=''%s'';\n', BatchID); % we pass the prefix; per-run jobid held inside .mat if needed
    fprintf(fid, 'run(''%s'');\n', post_process);
end

fclose(fid);

fprintf('\nMaster script created: %s\n', scriptName);
fprintf('Spawned %d parametric job scripts (%d originals × %d combos).\n', generatedCount, numel(files), numel(combKeys));

%% 8. EXECUTE THE MASTER SCRIPT ----------------------------------------------------------------------
warning('on');
run(scriptName);
fprintf('\n\nBatch processing finished. Enjoy the combinatorial explosion like an erudite pyromaniac .\n');

%% Auxiliary functions --------------------------------------------------------------------------------
function [combKeys, suffixes, flatNames, flatSets] = expand_var_ranges(var_ranges)
% expand_var_ranges - converts "var_ranges" definition into a Cartesian-product list
% Inputs:
%   var_ranges : Kx2 cell array: {name, rhsStringOrCellOfStrings}
% Outputs:
%   combKeys : 1xN cell, each item is Mx2 cell {name, rhsString} for one combo
%   suffixes : 1xN cell, each item is '_<v1>_<v2>_...' with '.' replaced by '_' (raw strings)
%   flatNames: 1xK variable names (order preserved)
%   flatSets : 1xK cell, each a 1xMi cellstr of candidate RHS strings
%
% Why this much ceremony? Because engineers don’t test, they suffer.

    % Validate shape because future-you will screw this up at 3am
    if ~iscell(var_ranges) || size(var_ranges,2) ~= 2
        error('var_ranges must be a Kx2 cell array. Congratulations on failing kindergarten grids.');
    end

    K = size(var_ranges,1);
    flatNames = var_ranges(:,1).';
    flatSets  = cell(1,K);

    for k = 1:K
        nameK = var_ranges{k,1};
        valsK = var_ranges{k,2};
        if ~(ischar(nameK) || isstring(nameK))
            error('Parameter name at row %d is not a string. Do better.', k);
        end
        if ischar(valsK) || isstring(valsK)
            flatSets{k} = {char(valsK)};
        elseif iscell(valsK)
            if ~all(cellfun(@(c)ischar(c)||isstring(c), valsK))
                error('Row %d values must be strings of MATLAB code.', k);
            end
            flatSets{k} = cellfun(@char, valsK, 'UniformOutput', false);
        else
            error('Row %d: values must be string or cellstr. Not whatever this is.', k);
        end
    end

    % Cartesian product over indices, then map back to strings
    sizes = cellfun(@numel, flatSets);
    if any(sizes == 0)
        error('One of the parameter ranges is empty. That’s not a range, that’s a void.');
    end

    % Build grids of indices using ndgrid like adults
    idxGrids = cell(1,K);
    [idxGrids{:}] = ndgrid(flatSets{:}); %#ok<CCAT> we don’t want to scalar-expand strings

    % Unfortunately ndgrid over cells returns cells of cells. Flatten like a caveman.
    N = numel(idxGrids{1});
    combKeys = cell(1,N);
    suffixes = cell(1,N);

    for n = 1:N
        pairs = cell(0,2);
        suffix = '';
        for k = 1:K
            % Extract the kth choice string
            valk = idxGrids{k}{n};
            namek = char(flatNames{k});
            pairs(end+1,:) = {namek, valk}; %#ok<AGROW>
            suffix = [suffix, '_', sanitize_token(valk)]; %#ok<AGROW>
        end
        combKeys{n} = pairs;
        suffixes{n} = suffix;
    end
end

function token = sanitize_token(rhsString)
% sanitize_token - turn an RHS string of code into a raw suffix token
% Rules (per your royal decree):
%   - Print raw values only (no key names)
%   - Replace '.' with '_' because filesystems cry
%   - Keep everything else; you gave us strings, we echo your sins
% Edge cases:
%   - If rhsString is a literal string like '''abc''', we strip outer quotes to use abc
%   - If it’s a number-like token '1e-3', we use it as-is then '.'→'_'
%   - If it’s a function call, you’ll get something ugly. That’s the deal you wanted.
    s = char(rhsString);

    % If s starts and ends with single quotes (literal), strip them
    if numel(s) >= 2 && s(1)=='''' && s(end)==''''
        % collapse doubled quotes '' -> '
        inner = s(2:end-1);
        inner = regexprep(inner, '''''', '''');
        s = inner;
    end

    % Replace '.' with '_'
    s = strrep(s, '.', '_');

    % Also purge whitespace because nothing good comes from spaces in IDs
    s = regexprep(s, '\s+', '');

    % And maybe commas (if someone passes "1,2"), convert to '_' for determinism
    s = strrep(s, ',', '_');

    % Brackets out ([]) — you wanted raw, not math
    s = regexprep(s, '[\[\]\(\)]', '');

    % Finally, assassinate filesystem hostiles: slashes and colons
    s = regexprep(s, '[/\\:;]', '_');

    token = s;
end

function copy_job_script(srcPath, dstPath)
% copy_job_script - copies original job script to new per-combo script
    if ~isfile(srcPath)
        error('Source job script not found: %s', srcPath);
    end
    [dstDir,~,~] = fileparts(dstPath);
    if ~exist(dstDir,'dir')
        mkdir(dstDir);
    end
    [status,msg] = copyfile(srcPath, dstPath, 'f');
    if ~status
        error('Failed to copy job script to %s: %s', dstPath, msg);
    end
end

function override_params(filename, newvars)
% OVERRIDE_PARAMS - Override simulation parameters in MATLAB scripts
%
% Syntax:
%   override_params(filename, newvars)
%
% Inputs:
%   filename - String with the path to the .m file to modify
%   newvars  - Cell array with N rows x 2 cols where:
%              Col 1: Parameter name (string/char)
%              Col 2: New value (string/char)
%
% Example:
%   newvars = {
%       'runLineParametersCalc', 'logical(1)';
%       'runSimuTool', 'logical(0)';
%       'decmp_flag', '1'
%   };
%   override_params('my_simulation.m', newvars)
%
% Notes:
%   This function replaces entire lines that start with the parameter name
%   followed by whitespace and an equals sign or whitespace and equals sign.
%   It makes a backup of the original file with '.bak' extension.

% Input validation
if ~ischar(filename) && ~isstring(filename)
    error('Filename must be a string or character array');
end

if ~iscell(newvars) || size(newvars, 2) ~= 2
    error('newvars must be a cell array with exactly 2 columns');
end

try
    % Check if file exists
    if ~exist(filename, 'file')
        error('File not found: %s', filename);
    end

    % Read the file content
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Cannot open file: %s', filename);
    end

    content = cell(1000, 1); % Pre-allocate memory for file content
    lineCount = 0;

    % Read file line by line
    line = fgetl(fileID);
    while ischar(line)
        lineCount = lineCount + 1;
        if lineCount > length(content)
            content = [content; cell(1000, 1)]; % Expand preallocated space if needed
        end
        content{lineCount} = line;
        line = fgetl(fileID);
    end

    % Trim extra preallocated cells
    content = content(1:lineCount);

    fclose(fileID);

    % Create backup file
    backup_filename = [filename, '.bak'];
    fileID = fopen(backup_filename, 'w');
    if fileID == -1
        error('Cannot create backup file: %s', backup_filename);
    end

    for i = 1:lineCount
        fprintf(fileID, '%s\n', content{i});
    end

    fclose(fileID);

    % Process each parameter
    replacementCount = zeros(size(newvars, 1), 1);

    for i = 1:size(newvars, 1)
        paramName = char(newvars{i, 1});
        newValue = char(newvars{i, 2});

        % Create pattern to match lines like "paramName = value" or "paramName= value"
        pattern1 = ['^', paramName, '\s*='];
        pattern2 = ['^', paramName, '='];

        for j = 1:lineCount
            line = content{j};

            % Check if line starts with parameter name followed by equals sign
            if ~isempty(regexp(line, pattern1, 'once')) || ~isempty(regexp(line, pattern2, 'once'))
                % Replace line with new parameter definition
                content{j} = [paramName, ' = ', newValue, ';'];
                replacementCount(i) = replacementCount(i) + 1;
            end
        end
    end

    % Write modified content back to the file
    fileID = fopen(filename, 'w');
    if fileID == -1
        error('Cannot write to file: %s', filename);
    end

    for i = 1:lineCount
        fprintf(fileID, '%s\n', content{i});
    end

    fclose(fileID);

    % Report results
    disp(['File modified: ', filename]);
    disp(['Backup created: ', backup_filename]);

    % Report parameter replacements
    for i = 1:size(newvars, 1)
        paramName = char(newvars{i, 1});
        if replacementCount(i) == 0
            warning('Parameter "%s" was not found in the file', paramName);
        elseif replacementCount(i) > 1
            warning('Parameter "%s" was found and replaced %d times', paramName, replacementCount(i));
        else
            disp(['Parameter replaced: ', paramName]);
        end
    end

catch ME
    % Handle errors
    if exist('fileID', 'var') && fileID ~= -1
        fclose(fileID);
    end
    rethrow(ME);
end
end

function userString = getStringInput()
% Pre-allocate empty string for error cases
userString = '';

try
    % Create input dialog with error checking
    prompt = {'Enter file filter (default *.m):'};
    dlgTitle = 'File filter';
    defaultInput = {'*.m'};
    dims = [1 50]; % Dialog box dimensions [rows columns]

    % Display the input dialog and get result
    answer = inputdlg(prompt, dlgTitle, dims, defaultInput);

    % Check if user clicked Cancel or closed dialog
    if isempty(answer)
        warning('User cancelled the input dialog');
        return;
    end

    % Extract the string from the cell array
    userString = answer{1};

    % Validate input if needed
    if isempty(userString)
        warning('Empty string entered');
    end

catch ME
    warning('Error getting string input: %s', ME.message);
end
end

function filePath = getFilePath(initialDir, fileTypes)
% GETFILEPATH - Get file path using a dialog box
%
% Syntax:
%   filePath = getFilePath()
%   filePath = getFilePath(initialDir)
%   filePath = getFilePath(initialDir, fileTypes)
%
% Inputs:
%   initialDir - Optional initial directory (defaults to pwd)
%   fileTypes  - Optional cell array of file filter specs
%               e.g., {'*.m', 'MATLAB Files (*.m)'; '*.mat', 'MAT-Files (*.mat)'}
%
% Outputs:
%   filePath - Selected file path (empty if cancelled)

% Set defaults if not provided
if ~exist('initialDir', 'var') || isempty(initialDir)
    initialDir = pwd;
end

if ~exist('fileTypes', 'var') || isempty(fileTypes)
    fileTypes = {'*.m;*.mat', 'MATLAB Files (*.m, *.mat)';
        '*.txt', 'Text Files (*.txt)';
        '*.*', 'All Files (*.*)'};
end

% Pre-allocate empty result for error cases
filePath = '';

try
    % Show file selection dialog
    [fileName, pathName] = uigetfile(fileTypes, 'Select a file', initialDir);

    % Check if user clicked Cancel
    if isequal(fileName, 0) || isequal(pathName, 0)
        warning('File selection was cancelled');
        return;
    end

    % Combine path and filename
    filePath = fullfile(pathName, fileName);

    % Additional validation if needed
    if ~exist(filePath, 'file')
        warning('Selected file does not exist: %s', filePath);
        filePath = '';
    end

catch ME
    warning('Error in file selection: %s', ME.message);
end
end

function cleanup_directory(tgtdir)
    % Delete all .mat files in the target directory
    matFiles = fullfile(tgtdir, '*.mat');
    delete(matFiles);
    
    % Check if 'plots' folder exists and delete it if it does
    plotsDir = fullfile(tgtdir, 'plots');
    if exist(plotsDir, 'dir')
        rmdir(plotsDir, 's');  % 's' flag allows recursive deletion
    end
    
    % Error handling
    if ~isempty(dir(matFiles))
        warning('Some .mat files could not be deleted');
    end
    
    if exist(plotsDir, 'dir')
        error('Failed to delete plots directory');
    end
end

function [status, cmdout] = mapDrive(folderPath, driveLetter)
%MAPDRIVE Substitutes a drive letter for a folder path using the subst command.
%
%   SYNTAX:
%   [status, cmdout] = mapDrive(folderPath, driveLetter)
%
%   DESCRIPTION:
%   This function executes the Windows 'subst' command to associate a drive
%   letter with a specified folder path, creating a virtual drive. It
%   includes error checking and provides feedback on the command's success.
%
%   INPUTS:
%   folderPath  - A string containing the full path to the folder.
%                 (e.g., 'C:\Users\YourUser\Documents\MyProject')
%   driveLetter - A single character string for the desired drive letter.
%                 (e.g., 'Z')
%
%   OUTPUTS:
%   status      - The exit status of the system command. 0 indicates success.
%   cmdout      - The command output as a string, which can be useful for
%                 debugging if the command fails.
%
%   EXAMPLE:
%   % Map drive X: to a specific project folder
%   projectFolder = 'D:\Work\CurrentProject';
%   [status, ~] = mapDrive(projectFolder, 'X');
%   if status == 0
%       disp('Drive X: mapped successfully.');
%   end
%
%   % To unmap the drive later:
%   % system('subst X: /D');

% --- Input Validation ---
% Check if the correct number of arguments are provided.
if nargin ~= 2
    error('mapDrive:IncorrectInputCount', 'This function requires exactly two input arguments: folderPath and driveLetter.');
end
% Check if inputs are strings.
if ~ischar(folderPath) || ~ischar(driveLetter)
    error('mapDrive:InvalidInputType', 'Both folderPath and driveLetter must be character strings.');
end
% Check if the folder exists.
if ~isfolder(folderPath)
    error('mapDrive:FolderNotFound', 'The specified folder does not exist: %s', folderPath);
end


% --- Command Execution ---

% 1. Construct the full command string for the command prompt.
%    Quotation marks are added around the path to handle any spaces.
substCommand = sprintf('subst %s: "%s"', driveLetter, folderPath);

% (Optional) Display the command that will be executed.
fprintf('Executing command: %s\n', substCommand);

% 2. Execute the command using the system() function.
[status, cmdout] = system(substCommand);

% 3. Check the result and provide feedback to the user.
if status == 0
    fprintf('Success! Drive %s: is now mapped to "%s".\n', driveLetter, folderPath);
    fprintf('\nTo remove this virtual drive later, run:\n');
    fprintf('  system(''subst %s: /D'');\n', driveLetter);
else
    fprintf('Error: The command failed to execute.\n');
    fprintf('Details:\n%s\n', cmdout);
end

end
