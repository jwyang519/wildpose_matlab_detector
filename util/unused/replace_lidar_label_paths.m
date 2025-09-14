function replace_lidar_label_paths(label_folder, varargin)
% REPLACE_LIDAR_LABEL_PATHS - Replace specific path prefixes in LiDAR label files
%
% Usage:
%   replace_lidar_label_paths(label_folder)
%   replace_lidar_label_paths(label_folder, 'OldPaths', old_paths, 'NewPath', new_path)
%   replace_lidar_label_paths(label_folder, 'DryRun', true)
%   replace_lidar_label_paths(label_folder, 'FieldNames', {'DataSource', 'FilePath'})
%
% Parameters:
%   label_folder  - Folder containing .mat label files
%   OldPaths      - Cell array of old path prefixes to replace
%   NewPath       - New path prefix to use
%   FieldNames    - Cell array of field names to update (default: auto-detect)
%   DryRun        - If true, only show what would be changed (default: false)
%   BackupFolder  - Folder to save backups (default: create backup subfolder)
%   Recursive     - Search subfolders recursively (default: false)
%   NormalizePaths - Convert all backslashes to forward slashes (default: true)

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'label_folder', @(x) ischar(x) || isstring(x));
    
    % Default old paths based on your requirements
    default_old_paths = {'/Volumes/T7_Shield/v1.1_pcd_s/', ...
                        'D:\v1.1_pcd_s\', ...
                        'D:\WildPosev1.0\Cheetah\2022-06-13\'};
    default_new_path = '/Volumes/T7_Shield/WildPose/WildPose_Proj/data';
    
    addParameter(p, 'OldPaths', default_old_paths, @iscell);
    addParameter(p, 'NewPath', default_new_path, @(x) ischar(x) || isstring(x));
    addParameter(p, 'FieldNames', {}, @iscell);
    addParameter(p, 'DryRun', false, @islogical);
    addParameter(p, 'BackupFolder', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Recursive', false, @islogical);
    addParameter(p, 'NormalizePaths', true, @islogical);
    
    parse(p, label_folder, varargin{:});
    
    % Extract parameters
    old_paths = p.Results.OldPaths;
    new_path = p.Results.NewPath;
    field_names = p.Results.FieldNames;
    dry_run = p.Results.DryRun;
    backup_folder = p.Results.BackupFolder;
    recursive = p.Results.Recursive;
    normalize_paths = p.Results.NormalizePaths;
    
    % Normalize the new path
    if normalize_paths
        new_path = strrep(new_path, '\', '/');
        % Ensure new path ends with /
        if ~endsWith(new_path, '/')
            new_path = [new_path '/'];
        end
    end
    
    % Normalize old paths for consistent matching
    normalized_old_paths = cell(size(old_paths));
    for i = 1:length(old_paths)
        normalized_old_paths{i} = strrep(old_paths{i}, '\', '/');
        % Ensure old paths end with /
        if ~endsWith(normalized_old_paths{i}, '/')
            normalized_old_paths{i} = [normalized_old_paths{i} '/'];
        end
    end
    
    fprintf('=== PATH REPLACEMENT CONFIGURATION ===\n');
    fprintf('Old paths to replace:\n');
    for i = 1:length(old_paths)
        fprintf('  %d. "%s"\n', i, old_paths{i});
    end
    fprintf('New path: "%s"\n', new_path);
    fprintf('Normalize paths (\ to /): %s\n', mat2str(normalize_paths));
    fprintf('\n');
    
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
            fprintf('No path fields detected. Searching all string fields...\n');
            field_names = {'*'}; % Will search all fields
        else
            fprintf('Detected path fields: %s\n', strjoin(field_names, ', '));
        end
    end
    
    % Process each file
    success_count = 0;
    error_count = 0;
    total_replacements = 0;
    
    for i = 1:length(mat_files)
        file_path = fullfile(mat_files(i).folder, mat_files(i).name);
        
        try
            fprintf('Processing (%d/%d): %s\n', i, length(mat_files), mat_files(i).name);
            
            % Load the file
            data = load(file_path);
            
            % Track if any changes were made
            changes_made = false;
            change_log = {};
            file_replacements = 0;
            
            % Update paths in all variables
            var_names = fieldnames(data);
            for j = 1:length(var_names)
                [data.(var_names{j}), field_changes, field_log, field_count] = ...
                    replace_paths_in_variable(data.(var_names{j}), field_names, ...
                    old_paths, normalized_old_paths, new_path, normalize_paths);
                
                if field_changes
                    changes_made = true;
                    change_log = [change_log, field_log];
                    file_replacements = file_replacements + field_count;
                end
            end
            
            % Display changes
            if changes_made
                fprintf('  Changes made (%d replacements):\n', file_replacements);
                for k = 1:min(10, length(change_log)) % Show max 10 changes per file
                    fprintf('    %s\n', change_log{k});
                end
                if length(change_log) > 10
                    fprintf('    ... and %d more changes\n', length(change_log) - 10);
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
                total_replacements = total_replacements + file_replacements;
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
    fprintf('Total path replacements: %d\n', total_replacements);
    fprintf('Errors: %d\n', error_count);
    
    if dry_run
        fprintf('\nThis was a DRY RUN - no files were actually modified.\n');
        fprintf('Run with DryRun=false to apply changes.\n');
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
            
            if length(var) == 1
                if ischar(var.(current_name)) || isstring(var.(current_name))
                    str_val = char(var.(current_name));
                    if contains_path_indicators(str_val)
                        path_fields{end+1} = current_name;
                    end
                elseif isstruct(var.(current_name))
                    sub_fields = find_path_fields_recursive(var.(current_name), current_name);
                    path_fields = [path_fields, sub_fields];
                end
            end
        end
    end
end

function [updated_var, changes_made, change_log, replacement_count] = replace_paths_in_variable(var, field_names, old_paths, normalized_old_paths, new_path, normalize_paths)
    % Replace paths in a variable recursively
    
    updated_var = var;
    changes_made = false;
    change_log = {};
    replacement_count = 0;
    
    if isstruct(var)
        for i = 1:length(var)
            struct_fields = fieldnames(var);
            
            for j = 1:length(struct_fields)
                field_name = struct_fields{j};
                
                % Check if we should process this field
                should_process = false;
                if any(strcmp(field_names, '*'))
                    should_process = true;
                elseif any(strcmp(field_names, field_name))
                    should_process = true;
                end
                
                if should_process
                    old_value = var(i).(field_name);
                    
                    if ischar(old_value) || isstring(old_value)
                        [new_value, was_changed] = replace_path_string(char(old_value), old_paths, normalized_old_paths, new_path, normalize_paths);
                        
                        if was_changed
                            updated_var(i).(field_name) = new_value;
                            changes_made = true;
                            replacement_count = replacement_count + 1;
                            change_log{end+1} = sprintf('%s: "%s" -> "%s"', field_name, char(old_value), new_value);
                        end
                        
                    elseif iscell(old_value)
                        [new_cell, cell_changed, cell_log, cell_count] = replace_cell_paths(old_value, old_paths, normalized_old_paths, new_path, normalize_paths);
                        if cell_changed
                            updated_var(i).(field_name) = new_cell;
                            changes_made = true;
                            replacement_count = replacement_count + cell_count;
                            change_log = [change_log, cell_log];
                        end
                        
                    elseif isstruct(old_value)
                        % Recursively process nested structs
                        [nested_result, nested_changed, nested_log, nested_count] = replace_paths_in_variable(old_value, field_names, old_paths, normalized_old_paths, new_path, normalize_paths);
                        if nested_changed
                            updated_var(i).(field_name) = nested_result;
                            changes_made = true;
                            replacement_count = replacement_count + nested_count;
                            % Prefix nested logs with field name
                            for k = 1:length(nested_log)
                                change_log{end+1} = sprintf('%s.%s', field_name, nested_log{k});
                            end
                        end
                    end
                end
            end
        end
    end
end

function [updated_cell, changes_made, change_log, replacement_count] = replace_cell_paths(cell_var, old_paths, normalized_old_paths, new_path, normalize_paths)
    % Replace paths in cell arrays
    
    updated_cell = cell_var;
    changes_made = false;
    change_log = {};
    replacement_count = 0;
    
    for i = 1:numel(cell_var)
        if ischar(cell_var{i}) || isstring(cell_var{i})
            old_path = char(cell_var{i});
            [new_path_str, was_changed] = replace_path_string(old_path, old_paths, normalized_old_paths, new_path, normalize_paths);
            
            if was_changed
                updated_cell{i} = new_path_str;
                changes_made = true;
                replacement_count = replacement_count + 1;
                change_log{end+1} = sprintf('Cell %d: "%s" -> "%s"', i, old_path, new_path_str);
            end
        end
    end
end

function [new_string, was_changed] = replace_path_string(input_string, old_paths, normalized_old_paths, new_path, normalize_paths)
    % Replace path prefixes in a string
    
    new_string = input_string;
    was_changed = false;
    
    if isempty(input_string)
        return;
    end
    
    % First normalize the input string for comparison
    normalized_input = strrep(input_string, '\', '/');
    
    % Try to replace each old path prefix
    for i = 1:length(old_paths)
        old_path_original = old_paths{i};
        old_path_normalized = normalized_old_paths{i};
        
        % Check if the normalized input starts with this old path
        if startsWith(normalized_input, old_path_normalized, 'IgnoreCase', true)
            % Replace the prefix
            remaining_path = normalized_input(length(old_path_normalized) + 1:end);
            new_string = [new_path remaining_path];
            was_changed = true;
            break; % Only replace the first match
            
        elseif startsWith(input_string, old_path_original, 'IgnoreCase', true)
            % Direct match with original format
            remaining_path = input_string(length(old_path_original) + 1:end);
            if normalize_paths
                remaining_path = strrep(remaining_path, '\', '/');
            end
            new_string = [new_path remaining_path];
            was_changed = true;
            break;
        end
    end
    
    % If no prefix replacement but we want to normalize paths
    if ~was_changed && normalize_paths && contains(input_string, '\')
        new_string = strrep(input_string, '\', '/');
        was_changed = true;
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
                      'C:', 'D:', '/home/', '/data/', 'file://', 'http://', '/Volumes/'};
    
    for i = 1:length(path_indicators)
        if contains(str, path_indicators{i})
            has_path = true;
            break;
        end
    end
end