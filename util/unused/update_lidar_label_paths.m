function update_lidar_label_paths(label_folder, varargin)
% UPDATE_LIDAR_LABEL_PATHS - Update absolute paths to relative paths in LiDAR label files
%
% Usage:
%   update_lidar_label_paths(label_folder)
%   update_lidar_label_paths(label_folder, 'FieldNames', {'DataSource', 'FilePath'})
%   update_lidar_label_paths(label_folder, 'RootFolder', '/path/to/project/root')
%   update_lidar_label_paths(label_folder, 'DryRun', true)
%   update_lidar_label_paths(label_folder, 'BackupFolder', '/path/to/backup')
%
% Parameters:
%   label_folder  - Folder containing .mat label files
%   FieldNames    - Cell array of field names to update (default: auto-detect)
%   RootFolder    - Root folder for relative paths (default: label_folder)
%   DryRun        - If true, only show what would be changed (default: false)
%   BackupFolder  - Folder to save backups (default: create backup subfolder)
%   Recursive     - Search subfolders recursively (default: false)

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'label_folder', @(x) ischar(x) || isstring(x));
    addParameter(p, 'FieldNames', {}, @iscell);
    addParameter(p, 'RootFolder', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'DryRun', false, @islogical);
    addParameter(p, 'BackupFolder', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Recursive', false, @islogical);
    
    parse(p, label_folder, varargin{:});
    
    % Extract parameters
    field_names = p.Results.FieldNames;
    root_folder = p.Results.RootFolder;
    dry_run = p.Results.DryRun;
    backup_folder = p.Results.BackupFolder;
    recursive = p.Results.Recursive;
    
    % Set default root folder
    if isempty(root_folder)
        root_folder = label_folder;
    end
    
    % Find all .mat files
    if recursive
        mat_files = dir(fullfile(label_folder, '**', '*.mat'));
    else
        mat_files = dir(fullfile(label_folder, '*.mat'));
    end
    
    if isempty(mat_files)
        fprintf('No .mat files found in the specified folder.\n');
        return;
    end
    
    fprintf('Found %d .mat files to process\n', length(mat_files));
    
    % Create backup folder if not dry run
    if ~dry_run
        if isempty(backup_folder)
            backup_folder = fullfile(label_folder, 'backup_original');
        end
        if ~exist(backup_folder, 'dir')
            mkdir(backup_folder);
            fprintf('Created backup folder: %s\n', backup_folder);
        end
    end
    
    % Auto-detect field names if not provided
    if isempty(field_names)
        fprintf('Auto-detecting path fields...\n');
        field_names = auto_detect_path_fields(mat_files(1:min(3, end)), label_folder);
        if isempty(field_names)
            fprintf('No path fields detected. Please specify field names manually.\n');
            return;
        end
        fprintf('Detected path fields: %s\n', strjoin(field_names, ', '));
    end
    
    % Process each file
    success_count = 0;
    error_count = 0;
    
    for i = 1:length(mat_files)
        file_path = fullfile(mat_files(i).folder, mat_files(i).name);
        
        try
            fprintf('Processing (%d/%d): %s\n', i, length(mat_files), mat_files(i).name);
            
            % Load the file
            data = load(file_path);
            
            % Track if any changes were made
            changes_made = false;
            change_log = {};
            
            % Update paths in specified fields
            var_names = fieldnames(data);
            for j = 1:length(var_names)
                [data.(var_names{j}), field_changes, field_log] = ...
                    update_paths_in_variable(data.(var_names{j}), field_names, ...
                    root_folder, mat_files(i).folder);
                
                if field_changes
                    changes_made = true;
                    change_log = [change_log, field_log];
                end
            end
            
            % Display changes
            if changes_made
                fprintf('  Changes made:\n');
                for k = 1:length(change_log)
                    fprintf('    %s\n', change_log{k});
                end
                
                if ~dry_run
                    % Create backup
                    [~, name, ext] = fileparts(mat_files(i).name);
                    backup_path = fullfile(backup_folder, [name '_backup' ext]);
                    copyfile(file_path, backup_path);
                    
                    % Save updated file
                    save(file_path, '-struct', 'data');
                    fprintf('  File updated and backup saved\n');
                end
                success_count = success_count + 1;
            else
                fprintf('  No changes needed\n');
            end
            
        catch ME
            fprintf('  ERROR: %s\n', ME.message);
            error_count = error_count + 1;
        end
        
        fprintf('\n');
    end
    
    % Summary
    fprintf('=== PROCESSING COMPLETE ===\n');
    fprintf('Files processed: %d\n', length(mat_files));
    fprintf('Files updated: %d\n', success_count);
    fprintf('Errors: %d\n', error_count);
    
    if dry_run
        fprintf('\nThis was a DRY RUN - no files were actually modified.\n');
        fprintf('Set DryRun to false to apply changes.\n');
    end
end

function field_names = auto_detect_path_fields(mat_files, folder_path)
    % Auto-detect fields that contain file paths
    
    field_names = {};
    common_path_fields = {'DataSource', 'dataSource', 'data_source', ...
                         'FilePath', 'filepath', 'file_path', ...
                         'ImagePath', 'imagePath', 'image_path', ...
                         'PointCloudPath', 'pointCloudPath', 'pointcloud_path'};
    
    for i = 1:length(mat_files)
        file_path = fullfile(mat_files(i).folder, mat_files(i).name);
        
        try
            data = load(file_path);
            detected_fields = find_path_fields_recursive(data, '');
            field_names = union(field_names, detected_fields);
        catch
            continue;
        end
        
        % Also check for common field names
        var_names = fieldnames(data);
        for j = 1:length(var_names)
            if isstruct(data.(var_names{j}))
                struct_fields = fieldnames(data.(var_names{j}));
                common_found = intersect(struct_fields, common_path_fields);
                field_names = union(field_names, common_found);
            end
        end
    end
end

function path_fields = find_path_fields_recursive(var, prefix)
    % Recursively find fields containing paths
    
    path_fields = {};
    
    if isstruct(var)
        field_names = fieldnames(var);
        for i = 1:length(field_names)
            current_name = field_names{i};
            if ~isempty(prefix)
                full_name = [prefix '.' current_name];
            else
                full_name = current_name;
            end
            
            if length(var) == 1
                if ischar(var.(current_name)) || isstring(var.(current_name))
                    str_val = char(var.(current_name));
                    if contains_path_indicators(str_val)
                        path_fields{end+1} = current_name;
                    end
                elseif isstruct(var.(current_name))
                    sub_fields = find_path_fields_recursive(var.(current_name), full_name);
                    path_fields = [path_fields, sub_fields];
                end
            end
        end
    elseif iscell(var)
        for i = 1:min(5, length(var))
            if ischar(var{i}) || isstring(var{i})
                if contains_path_indicators(char(var{i}))
                    path_fields{end+1} = 'cell_content';
                    break;
                end
            end
        end
    end
end

function [updated_var, changes_made, change_log] = update_paths_in_variable(var, field_names, root_folder, current_folder)
    % Update paths in a variable recursively
    
    updated_var = var;
    changes_made = false;
    change_log = {};
    
    if isstruct(var)
        for i = 1:length(var)
            for j = 1:length(field_names)
                field_name = field_names{j};
                
                if isfield(var, field_name)
                    old_value = var(i).(field_name);
                    
                    if ischar(old_value) || isstring(old_value)
                        new_value = convert_to_relative_path(char(old_value), root_folder, current_folder);
                        
                        if ~strcmp(old_value, new_value)
                            updated_var(i).(field_name) = new_value;
                            changes_made = true;
                            change_log{end+1} = sprintf('%s: %s -> %s', field_name, char(old_value), new_value);
                        end
                    elseif iscell(old_value)
                        [new_cell, cell_changed, cell_log] = update_cell_paths(old_value, root_folder, current_folder);
                        if cell_changed
                            updated_var(i).(field_name) = new_cell;
                            changes_made = true;
                            change_log = [change_log, cell_log];
                        end
                    end
                end
            end
        end
    end
end

function [updated_cell, changes_made, change_log] = update_cell_paths(cell_var, root_folder, current_folder)
    % Update paths in cell arrays
    
    updated_cell = cell_var;
    changes_made = false;
    change_log = {};
    
    for i = 1:numel(cell_var)
        if ischar(cell_var{i}) || isstring(cell_var{i})
            old_path = char(cell_var{i});
            new_path = convert_to_relative_path(old_path, root_folder, current_folder);
            
            if ~strcmp(old_path, new_path)
                updated_cell{i} = new_path;
                changes_made = true;
                change_log{end+1} = sprintf('Cell %d: %s -> %s', i, old_path, new_path);
            end
        end
    end
end

function relative_path = convert_to_relative_path(absolute_path, root_folder, current_folder)
    % Convert absolute path to relative path
    
    if isempty(absolute_path) || ~contains_path_indicators(absolute_path)
        relative_path = absolute_path;
        return;
    end
    
    % Try to make path relative to root folder first
    try
        relative_path = get_relative_path(absolute_path, root_folder);
        
        % If the relative path goes up too many levels, try relative to current folder
        if count(relative_path, '../') > 3
            alt_relative = get_relative_path(absolute_path, current_folder);
            if count(alt_relative, '../') < count(relative_path, '../')
                relative_path = alt_relative;
            end
        end
    catch
        relative_path = absolute_path;  % Keep original if conversion fails
    end
end

function rel_path = get_relative_path(target_path, base_path)
    % Get relative path from base to target
    
    % Normalize paths
    target_path = strrep(target_path, '\', '/');
    base_path = strrep(base_path, '\', '/');
    
    % Split paths into components
    target_parts = strsplit(target_path, '/');
    base_parts = strsplit(base_path, '/');
    
    % Remove empty parts
    target_parts = target_parts(~cellfun(@isempty, target_parts));
    base_parts = base_parts(~cellfun(@isempty, base_parts));
    
    % Find common prefix
    common_length = 0;
    for i = 1:min(length(target_parts), length(base_parts))
        if strcmpi(target_parts{i}, base_parts{i})
            common_length = i;
        else
            break;
        end
    end
    
    % Build relative path
    up_levels = length(base_parts) - common_length;
    down_parts = target_parts(common_length + 1:end);
    
    rel_parts = [repmat({'../'}, 1, up_levels), down_parts];
    
    if isempty(rel_parts)
        rel_path = './';
    else
        rel_path = strjoin(rel_parts, '/');
        rel_path = strrep(rel_path, '//', '/');
    end
end

function has_path = contains_path_indicators(str)
    % Check if a string contains path-like indicators
    has_path = false;
    
    if isempty(str)
        return;
    end
    
    % Common path indicators
    path_indicators = {'\', '/', '.mat', '.pcd', '.bin', '.png', '.jpg', '.tif', ...
                      'C:', 'D:', '/home/', '/data/', 'file://', 'http://'};
    
    for i = 1:length(path_indicators)
        if contains(str, path_indicators{i})
            has_path = true;
            break;
        end
    end
end