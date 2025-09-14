% matlab_path_workaround.m
% Workaround for MATLAB's stubborn changeFilePaths function

function update_paths_workaround(labelDir)
    if nargin < 1
        labelDir = "/Volumes/T7_Shield/WildPose/WildPose_Proj/label/28-8-2025_30/bad";
    end
    
    matFiles = dir(fullfile(labelDir, "*.mat"));
    
    % Define path mappings
    oldPrefixes = {'D:\v1.1_pcd_s\', ...
                  '/Volumes/T7_Shield/v1.1_pcd_s/', ...
                  'D:\WildPosev1.0\Cheetah\2022-06-13\'};
    newBase = '/Volumes/T7_Shield/WildPose/WildPose_Proj/data/';
    
    fprintf('=== MATLAB PATH UPDATE WORKAROUND ===\n');
    fprintf('Processing %d .mat files...\n', numel(matFiles));
    
    successCount = 0;
    errorCount = 0;
    
    for k = 1:numel(matFiles)
        matFile = fullfile(matFiles(k).folder, matFiles(k).name);
        
        try
            fprintf('\nProcessing (%d/%d): %s\n', k, numel(matFiles), matFiles(k).name);
            
            % Method 1: Try the "proper" changeFilePaths approach with exact matching
            if try_change_file_paths_method(matFile, oldPrefixes, newBase)
                fprintf('  ✓ Success with changeFilePaths method\n');
                successCount = successCount + 1;
                continue;
            end
            
            % Method 2: Manual reconstruction approach
            if try_manual_reconstruction_method(matFile, oldPrefixes, newBase)
                fprintf('  ✓ Success with manual reconstruction method\n');
                successCount = successCount + 1;
                continue;
            end
            
            % Method 3: Raw data manipulation (last resort)
            if try_raw_data_method(matFile, oldPrefixes, newBase)
                fprintf('  ✓ Success with raw data manipulation method\n');
                successCount = successCount + 1;
                continue;
            end
            
            fprintf('  ✗ All methods failed\n');
            errorCount = errorCount + 1;
            
        catch ME
            fprintf('  ✗ ERROR: %s\n', ME.message);
            errorCount = errorCount + 1;
        end
    end
    
    fprintf('\n=== SUMMARY ===\n');
    fprintf('Total files: %d\n', numel(matFiles));
    fprintf('Successfully updated: %d\n', successCount);
    fprintf('Failed: %d\n', errorCount);
end

function success = try_change_file_paths_method(matFile, oldPrefixes, newBase)
    success = false;
    
    % Load the file
    S = load(matFile, 'gTruth');
    originalSourceName = S.gTruth.DataSource.SourceName;
    
    fprintf('  Method 1: Trying changeFilePaths with directory creation...\n');
    fprintf('    Original: "%s"\n', originalSourceName);
    
    % Find matching prefix
    matchingPrefix = '';
    newPath = '';
    
    for i = 1:length(oldPrefixes)
        if startsWith(char(originalSourceName), oldPrefixes{i}, 'IgnoreCase', true)
            matchingPrefix = oldPrefixes{i};
            remainingPath = char(originalSourceName(length(oldPrefixes{i})+1:end));
            remainingPath = strrep(remainingPath, '\', '/');
            newPath = [newBase, remainingPath];
            break;
        end
    end
    
    if isempty(matchingPrefix)
        fprintf('    No matching prefix found\n');
        return;
    end
    
    fprintf('    Expected new path: "%s"\n', newPath);
    
    % Create the directory structure if it doesn't exist
    newDir = fileparts(newPath);
    if ~exist(newDir, 'dir')
        fprintf('    Creating directory: "%s"\n', newDir);
        try
            mkdir(newDir);
        catch
            fprintf('    Failed to create directory\n');
            return;
        end
    end
    
    % Try changeFilePaths with multiple formats
    formats_to_try = {
        {matchingPrefix, newPath},                           % Exact path mapping
        {char(originalSourceName), newPath},                 % Full path mapping
        {matchingPrefix, newBase}                            % Prefix mapping
    };
    
    for fmt = 1:length(formats_to_try)
        altPaths = formats_to_try{fmt};
        fprintf('    Trying format %d: {"%s", "%s"}\n', fmt, altPaths{1}, altPaths{2});
        
        % Make a copy to test
        tempS = S;
        try
            unresolved = changeFilePaths(tempS.gTruth, altPaths);
            newSourceName = tempS.gTruth.DataSource.SourceName;
            
            if ~strcmp(char(originalSourceName), char(newSourceName))
                fprintf('    SUCCESS! Path changed to: "%s"\n', newSourceName);
                % Save the successful change
                save(matFile, '-struct', 'tempS');
                success = true;
                return;
            else
                fprintf('    No change in path\n');
            end
        catch ME
            fprintf('    Format %d failed: %s\n', fmt, ME.message);
        end
    end
end

function success = try_manual_reconstruction_method(matFile, oldPrefixes, newBase)
    success = false;
    
    fprintf('  Method 2: Trying manual reconstruction...\n');
    
    try
        % Load all variables from the .mat file
        S = load(matFile);
        
        if ~isfield(S, 'gTruth')
            fprintf('    No gTruth field found\n');
            return;
        end
        
        originalSourceName = S.gTruth.DataSource.SourceName;
        fprintf('    Original path: "%s"\n', originalSourceName);
        
        % Find the matching prefix
        newPath = '';
        for i = 1:length(oldPrefixes)
            if startsWith(char(originalSourceName), oldPrefixes{i}, 'IgnoreCase', true)
                remainingPath = char(originalSourceName(length(oldPrefixes{i})+1:end));
                remainingPath = strrep(remainingPath, '\', '/');
                newPath = [newBase, remainingPath];
                break;
            end
        end
        
        if isempty(newPath)
            fprintf('    No matching prefix found\n');
            return;
        end
        
        fprintf('    Target path: "%s"\n', newPath);
        
        % Create a new groundTruthLidar object with the updated path
        % We'll try to construct a new DataSource with the correct path
        
        % Get the DataSource class and try to create a new one
        oldDS = S.gTruth.DataSource;
        dsClass = class(oldDS);
        
        fprintf('    DataSource class: %s\n', dsClass);
        
        % Try to create new DataSource with updated path
        try
            switch dsClass
                case 'vision.labeler.loading.PointCloudSequenceSource'
                    % Create new PointCloudSequenceSource with new path
                    newDS = vision.labeler.loading.PointCloudSequenceSource(newPath);
                    
                case 'vision.labeler.loading.VelodyneLidarSource' 
                    newDS = vision.labeler.loading.VelodyneLidarSource(newPath);
                    
                otherwise
                    fprintf('    Unsupported DataSource class: %s\n', dsClass);
                    return;
            end
            
            % Create new groundTruthLidar with updated DataSource
            newGTruth = groundTruthLidar(newDS, S.gTruth.LabelDefinitions, S.gTruth.LabelData);
            
            % Replace the gTruth in the structure
            S.gTruth = newGTruth;
            
            % Save the updated file
            save(matFile, '-struct', 'S');
            
            fprintf('    Successfully reconstructed with new path\n');
            success = true;
            
        catch ME
            fprintf('    Reconstruction failed: %s\n', ME.message);
        end
        
    catch ME
        fprintf('    Method 2 error: %s\n', ME.message);
    end
end

function success = try_raw_data_method(matFile, oldPrefixes, newBase)
    success = false;
    
    fprintf('  Method 3: Trying raw data manipulation...\n');
    
    try
        % Read the .mat file as raw data
        data = load(matFile);
        
        if ~isfield(data, 'gTruth')
            return;
        end
        
        % Convert to struct to access private properties
        gTruthStruct = struct(data.gTruth);
        
        fprintf('    Available fields: %s\n', strjoin(fieldnames(gTruthStruct), ', '));
        
        % Look for the private DSource property
        if isfield(gTruthStruct, 'DSource')
            dsStruct = struct(gTruthStruct.DSource);
            fprintf('    DSource fields: %s\n', strjoin(fieldnames(dsStruct), ', '));
            
            % Try to modify the SourceName in the struct
            originalPath = dsStruct.SourceName;
            fprintf('    Original DSource path: "%s"\n', originalPath);
            
            % Find and replace the path
            for i = 1:length(oldPrefixes)
                if startsWith(char(originalPath), oldPrefixes{i}, 'IgnoreCase', true)
                    remainingPath = char(originalPath(length(oldPrefixes{i})+1:end));
                    remainingPath = strrep(remainingPath, '\', '/');
                    newPath = [newBase, remainingPath];
                    
                    fprintf('    Attempting to set new path: "%s"\n', newPath);
                    
                    % This is a desperate attempt - probably won't work
                    % but worth trying
                    try
                        dsStruct.SourceName = newPath;
                        dsStruct.IsLoaded = false;
                        
                        % Try to save back (this will likely fail)
                        gTruthStruct.DSource = dsStruct;
                        data.gTruth = gTruthStruct;
                        
                        save(matFile, '-struct', 'data');
                        
                        fprintf('    Raw manipulation succeeded!\n');
                        success = true;
                        return;
                        
                    catch ME
                        fprintf('    Raw manipulation failed: %s\n', ME.message);
                    end
                    break;
                end
            end
        end
        
    catch ME
        fprintf('    Method 3 error: %s\n', ME.message);
    end
end