function inspect_lidar_labels(folder_path)
% INSPECT_LIDAR_LABELS - Dry run to inspect LiDAR label .mat files
% This function examines the structure of your .mat files to understand
% what fields contain paths that need to be updated
%
% Usage: inspect_lidar_labels('path/to/your/label/folder')

    if nargin < 1
        folder_path = uigetdir(pwd, 'Select folder containing .mat label files');
        if folder_path == 0
            return;
        end
    end
    
    % Find all .mat files in the folder
    mat_files = dir(fullfile(folder_path, '*.mat'));
    
    if isempty(mat_files)
        fprintf('No .mat files found in the specified folder.\n');
        return;
    end
    
    fprintf('Found %d .mat files to inspect\n\n', length(mat_files));
    
    % Inspect first few files to understand structure
    num_to_inspect = min(10, length(mat_files));
    
    for i = 1:num_to_inspect
        file_path = fullfile(folder_path, mat_files(i).name);
        fprintf('=== Inspecting file %d: %s ===\n', i, mat_files(i).name);
        
        try
            % Load the .mat file
            data = load(file_path);
            
            % Display top-level variables
            var_names = fieldnames(data);
            fprintf('Top-level variables: %s\n', strjoin(var_names, ', '));
            
            % Recursively inspect each variable for potential path fields
            for j = 1:length(var_names)
                fprintf('\n--- Variable: %s ---\n', var_names{j});
                inspect_variable(data.(var_names{j}), var_names{j}, 0);
            end
            
        catch ME
            fprintf('Error loading file: %s\n', ME.message);
        end
        
        fprintf('\n');
    end
    
    % Summary
    fprintf('=== INSPECTION SUMMARY ===\n');
    fprintf('Common path-related fields to look for:\n');
    fprintf('- DataSource, dataSource, data_source\n');
    fprintf('- FilePath, filepath, file_path\n');
    fprintf('- ImagePath, imagePath, image_path\n');
    fprintf('- PointCloudPath, pointCloudPath\n');
    fprintf('- Any field containing file paths or URLs\n\n');
    
    fprintf('Next steps:\n');
    fprintf('1. Review the output above to identify which fields contain paths\n');
    fprintf('2. Use update_lidar_label_paths() function with appropriate field names\n');
    fprintf('3. Test on a backup copy first!\n');
end

function inspect_variable(var, var_name, depth)
    % Recursive function to inspect variables for potential path fields
    
    indent = repmat('  ', 1, depth);
    
    if isstruct(var)
        if length(var) > 1
            fprintf('%s%s: struct array [%s]\n', indent, var_name, mat2str(size(var)));
            if length(var) <= 5  % Only inspect first few elements
                for k = 1:min(length(var), 2)
                    fprintf('%s  Element %d:\n', indent, k);
                    field_names = fieldnames(var(k));
                    for f = 1:length(field_names)
                        inspect_variable(var(k).(field_names{f}), field_names{f}, depth + 2);
                    end
                end
            end
        else
            fprintf('%s%s: struct\n', indent, var_name);
            field_names = fieldnames(var);
            for f = 1:length(field_names)
                inspect_variable(var.(field_names{f}), field_names{f}, depth + 1);
            end
        end
        
    elseif iscell(var)
        fprintf('%s%s: cell array [%s]\n', indent, var_name, mat2str(size(var)));
        % Check if cells contain strings that might be paths
        if ~isempty(var)
            sample_cells = var(1:min(3, numel(var)));
            for c = 1:length(sample_cells)
                if ischar(sample_cells{c}) || isstring(sample_cells{c})
                    str_val = char(sample_cells{c});
                    if contains_path_indicators(str_val)
                        fprintf('%s  Cell %d (POTENTIAL PATH): %s\n', indent, c, str_val);
                    elseif length(str_val) < 100  % Only show short strings
                        fprintf('%s  Cell %d: %s\n', indent, c, str_val);
                    else
                        fprintf('%s  Cell %d: <long string>\n', indent, c);
                    end
                end
            end
        end
        
    elseif ischar(var) || isstring(var)
        str_val = char(var);
        if contains_path_indicators(str_val)
            fprintf('%s%s (POTENTIAL PATH): %s\n', indent, var_name, str_val);
        elseif length(str_val) < 100
            fprintf('%s%s: %s\n', indent, var_name, str_val);
        else
            fprintf('%s%s: <long string>\n', indent, var_name);
        end
        
    else
        % Numeric or other data types
        fprintf('%s%s: %s [%s]\n', indent, var_name, class(var), mat2str(size(var)));
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