% =========================================================================
% ParametricRunBuilder.m
%
% DESCRIPTION:
%   This script automates the process of running a series of parametric
%   simulations. It recursively finds "job" files (e.g., 'RunJob_*.m') in a
%   source directory, modifies specific parameters within them, and then
%   generates a master script ('RunAllJobs.m') to execute them sequentially.
%   After each job, it copies the resulting data and plots to a target
%   directory. Finally, it can run a post-processing script on all the
%   aggregated results.
%
% HOW TO USE:
%   1. Configure the parameters in the "CONFIGURATION" section below.
%   2. Set the source/target directories and the job file name pattern.
%   3. Run this script. It will automatically generate and execute
%      'RunAllJobs.m' in the current directory.
%
% =========================================================================

%% OPTIONAL - MAP FOLDER TO VIRTUAL DRIVE
% -------------------------------------------------------------------------
% Maps a folder with a long path to a virtual drive to make it shorter.
driveLetter = 'Z';
folderPath = '';    % Leave empty if not used

% Check if the folderPath variable is defined, not empty, and is a valid directory
if ~isempty(folderPath) && isfolder(folderPath)
    fprintf('Attempting to map drive %s: to path: %s\n', driveLetter, folderPath);
    % Call the function to map the drive
    [status, cmdout] = mapDrive(folderPath, driveLetter);
else
    fprintf('Drive mapping skipped: folderPath is empty or not a valid folder.\n');
end

%% 1. INITIAL WORKSPACE CLEANUP
% -------------------------------------------------------------------------
% Clears the workspace to ensure a clean run without old data interfering.
close all;           % Close all open figures
clear all;           % Clear all variables from the workspace
delete(allchild(0)); % Delete hidden graphics objects
clc;                 % Clear the command window

%% 2. CONFIGURATION
% -------------------------------------------------------------------------
% This section contains all the user-definable settings for the batch run.

% A unique identifier for this batch of runs. Used for post-processing.
BatchID = 'MyParametricRun';

% --- Parameters to Override in Each Job File ---
% This cell array defines the variables that will be changed in each of the
% job files found by the script.
% Format: {'VariableName', 'NewValueAsString'; ...}
newvars = {
    'runLineParametersCalc', 'logical(1)';
    'runSimuTool', 'logical(1)';
    'ZYprnt', 'logical(1)';
    'Modprnt', 'logical(1)';
    'f', 'transpose(logspaceppd(-1,6,20))';
    'ZYsave', 'logical(1)';
    'Figsave', 'logical(1)';
    'export2EMTP', 'logical(0)';
    'export2XML', 'logical(0)';
    'export2ULM', 'logical(0)';
    'decmp_flag', '9';
    };

% --- Directory and File Paths ---
% NOTE: If these variables are commented out or not defined, the script
% will prompt you to select them using a file dialog.

% PRO-TIP: Use the virtual drive Z:\ (or your preferred letter) to keep
% the scripts short.

% The root directory where the script will start searching for job files.
% srcdir = '/path/to/source/files';

% The directory where all results (.mat files, plots) will be copied.
% tgtdir = '/path/to/results/storage';

% The file name pattern for the job files to find (wildcards accepted).
% fname = 'RunJob_*1000*.m';

% (Optional) Path to a script that will run after all jobs are complete.
% post_process = '/path/MyNicePlots.m';

% --- Cleanup Flag ---
% If set to true, the script will delete the contents of the target
% directory ('tgtdir') and any job-specific subdirectories before starting.
% USE WITH CAUTION.
force_cleanup = true;


%% 3. SETUP & VALIDATION
% -------------------------------------------------------------------------
% This section prepares the environment, validates paths, and ensures
% all necessary information is available before proceeding.

% If srcdir is not defined, prompt the user to select it.
if ~exist('srcdir', 'var')
    srcdir = uigetdir('', 'Select source directory for recursive search');
    if srcdir == 0; error('No directory selected. Exiting...'); end
end

% Change the current directory to the source directory.
cd(srcdir);

% If tgtdir is not defined, prompt the user to select it.
if ~exist('tgtdir', 'var')
    tgtdir = uigetdir('', 'Select target directory for copied files');
    if tgtdir == 0; error('No target directory selected. Exiting...'); end
end

% If the post-processing script path is not defined, prompt the user.
if ~exist('post_process', 'var')
    post_process = getFilePath(); % Assumes getFilePath is a helper function
end

% If the job file name pattern is not defined, prompt the user.
if ~exist('fname', 'var')
    fname = getStringInput(); % Assumes getStringInput is a helper function
    if isempty(fname); fname = '*.m'; end
end

% Perform initial cleanup of the target directory if requested.
if force_cleanup
    % This assumes you have a helper function 'cleanup_directory'.
    cleanup_directory(tgtdir);
end


%% 4. FIND JOB FILES
% -------------------------------------------------------------------------
% Recursively search the source directory for files matching the pattern.

fprintf('Searching for files matching "%s" in "%s"...\n', fname, srcdir);
files = dir(fullfile(srcdir, '**', fname));

if isempty(files)
    warning('No matching files found in %s. Exiting.', srcdir);
    return;
end
fprintf('Found %d job files. Proceeding to generate execution script...\n', length(files));


%% 5. GENERATE THE MASTER EXECUTION SCRIPT
% -------------------------------------------------------------------------
% Create a new .m file ('RunAllJobs.m') that will contain the commands to
% run each job file in sequence.

scriptName = 'RunAllJobs.m';
fid = fopen(scriptName, 'w'); % Open the file for writing

% Write a header to the auto-generated script.
fprintf(fid, '%% Auto-generated script to run all jobs sequentially\n');
fprintf(fid, '%% Generated by ParametricRunBuilder.m on %s\n\n', datestr(now));

% Loop through each file found.
for i = 1:length(files)
    currentFile = fullfile(files(i).folder, files(i).name);

    % --- Pre-processing for each job ---
    % Modify the parameters in the original job file.
    % This assumes you have a helper function 'override_params'.
    override_params(currentFile, newvars);

    % Clean up the specific job's folder if requested.
    if force_cleanup
        % This assumes you have a helper function 'cleanup_directory'.
        cleanup_directory(files(i).folder);
    end

    % --- Write job execution commands to the script ---
    fprintf(fid, '%% --- Job %d/%d: %s ---\n', i, length(files), files(i).name);
    fprintf(fid, 'cd(''%s'');\n', files(i).folder); % Change to the job's directory
    fprintf(fid, 'try\n');
    fprintf(fid, '    fprintf(''\\n\\n======================================\\n'');\n');
    fprintf(fid, '    fprintf(''Now running job %d/%d: %s\\n'');\n', i, length(files), files(i).name);
    fprintf(fid, '    fprintf(''======================================\\n\\n'');\n');
    fprintf(fid, '    run(''%s'');\n', files(i).name); % Run the job script

    % --- Write result-copying commands to the script ---
    fprintf(fid, '    if exist(''jobid'', ''var'')\n');
    fprintf(fid, '        copyfile([jobid ''_NLTSimu.mat''], ''%s'', ''f'');\n', tgtdir);
    fprintf(fid, '        fprintf(''Copied %%s_NLTSimu.mat to target directory.\\n'', jobid);\n');
    fprintf(fid, '        copyfile([jobid ''_LineParam.mat''], ''%s'', ''f'');\n', tgtdir);
    fprintf(fid, '        fprintf(''Copied %%s_LineParam.mat to target directory.\\n'', jobid);\n');
    fprintf(fid, '        copyfile(''./plots'', fullfile(''%s'',''plots''));\n', tgtdir);
    fprintf(fid, '        fprintf(''Copied plots subfolder to target directory.\\n'');\n');
    fprintf(fid, '    else\n');
    fprintf(fid, '        warning(''jobid variable not defined after running %s. Results may not have been copied.'');\n', files(i).name);
    fprintf(fid, '    end\n');

    % --- Write error handling to the script ---
    fprintf(fid, 'catch ME\n');
    fprintf(fid, '    warning(''An error occurred in job %s: %%s'', ME.message);\n', files(i).name);
    fprintf(fid, 'end\n\n');
end

%% 6. ADD POST-PROCESSING STEP TO SCRIPT
% -------------------------------------------------------------------------
% If a post-processing script was defined, add it to the end of the master
% script.

fprintf(fid, 'fprintf(''\\n\\nAll jobs completed!\\n'');\n');

if ~isempty(post_process) && isfile(post_process)
    fprintf('Adding post-processing step...\n');
    fprintf(fid, 'fprintf(''Running post-processing script...\\n'');\n');
    fprintf(fid, 'batch_tag=''%s'';\n', BatchID); % Pass the batch tag
    fprintf(fid, 'run(''%s'');\n', post_process);
end

% Close the file handle for the master script.
fclose(fid);


%% 7. EXECUTE THE GENERATED SCRIPT
% -------------------------------------------------------------------------
% Finally, run the 'RunAllJobs.m' script that was just created.

fprintf('\nMaster script created: %s\n', scriptName);
fprintf('Executing script now to process all %d jobs...\n', length(files));

warning('on'); % Ensure warnings are enabled to catch potential issues
run(scriptName);

fprintf('\n\nBatch processing finished.\n');

%% Auxiliary functions --------------------------------------------------------------------------------
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
