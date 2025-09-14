% direct_path_update.m
% Directly modify the DataSource.SourceName property, bypassing changeFilePaths

labelDir = "/Volumes/T7_Shield/WildPose/WildPose_Proj/label/28-8-2025_30/bad";
matFiles = dir(fullfile(labelDir, "*.mat"));

% Define path mappings
oldPrefixes = {'D:\v1.1_pcd_s\', ...
              '/Volumes/T7_Shield/v1.1_pcd_s/', ...
              'D:\WildPosev1.0\Cheetah\2022-06-13\'};
newBase = '/Volumes/T7_Shield/WildPose/WildPose_Proj/data/';

fprintf('=== DIRECT PATH UPDATE ===\n');
fprintf('Processing %d .mat files...\n', numel(matFiles));
fprintf('Path mappings:\n');
for i = 1:length(oldPrefixes)
    fprintf('  "%s" -> "%s"\n', oldPrefixes{i}, newBase);
end
fprintf('\n');

successCount = 0;
errorCount = 0;
noChangeCount = 0;

for k = 1:numel(matFiles)
    matFile = fullfile(matFiles(k).folder, matFiles(k).name);
    
    try
        fprintf('Processing (%d/%d): %s\n', k, numel(matFiles), matFiles(k).name);
        
        % Load the file
        S = load(matFile, 'gTruth');
        gTruth = S.gTruth;
        
        % Track if any changes were made
        changesMade = false;
        
        % Process each DataSource
        for dsIdx = 1:length(gTruth.DataSource)
            originalPath = char(gTruth.DataSource(dsIdx).SourceName);
            fprintf('  DataSource %d: "%s"\n', dsIdx, originalPath);
            
            % Find matching prefix and replace
            newPath = originalPath;
            for prefixIdx = 1:length(oldPrefixes)
                oldPrefix = oldPrefixes{prefixIdx};
                
                if startsWith(originalPath, oldPrefix, 'IgnoreCase', true)
                    % Extract the remaining path after the prefix
                    remainingPath = originalPath(length(oldPrefix)+1:end);
                    
                    % Normalize path separators to forward slashes
                    remainingPath = strrep(remainingPath, '\', '/');
                    
                    % Construct new path
                    newPath = [newBase, remainingPath];
                    
                    fprintf('    -> Updating to: "%s"\n', newPath);
                    
                    % Directly modify the SourceName
                    gTruth.DataSource(dsIdx).SourceName = newPath;
                    
                    % Reset IsLoaded to false since path changed
                    gTruth.DataSource(dsIdx).IsLoaded = false;
                    
                    changesMade = true;
                    break;
                end
            end
            
            if strcmp(originalPath, newPath)
                fprintf('    -> No change needed\n');
            end
        end
        
        % Handle VoxelLabelData if present
        if istable(gTruth.LabelData) && any(strcmp(gTruth.LabelData.Properties.VariableNames, 'VoxelLabelData'))
            fprintf('  Processing VoxelLabelData...\n');
            voxelData = gTruth.LabelData.VoxelLabelData;
            
            for vIdx = 1:length(voxelData)
                if ~isempty(voxelData{vIdx}) && (ischar(voxelData{vIdx}) || isstring(voxelData{vIdx}))
                    originalVoxelPath = char(voxelData{vIdx});
                    
                    for prefixIdx = 1:length(oldPrefixes)
                        oldPrefix = oldPrefixes{prefixIdx};
                        
                        if startsWith(originalVoxelPath, oldPrefix, 'IgnoreCase', true)
                            remainingPath = originalVoxelPath(length(oldPrefix)+1:end);
                            remainingPath = strrep(remainingPath, '\', '/');
                            newVoxelPath = [newBase, remainingPath];
                            
                            fprintf('    VoxelData %d: "%s" -> "%s"\n', vIdx, originalVoxelPath, newVoxelPath);
                            
                            gTruth.LabelData.VoxelLabelData{vIdx} = newVoxelPath;
                            changesMade = true;
                            break;
                        end
                    end
                end
            end
        end
        
        if changesMade
            % Save the updated file
            save(matFile, '-struct', 'S', 'gTruth');
            fprintf('  ✓ Successfully updated and saved\n');
            successCount = successCount + 1;
        else
            fprintf('  - No changes made\n');
            noChangeCount = noChangeCount + 1;
        end
        
    catch ME
        fprintf('  ✗ ERROR: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('    at %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        errorCount = errorCount + 1;
    end
    
    fprintf('\n');
end

fprintf('=== SUMMARY ===\n');
fprintf('Total files processed: %d\n', numel(matFiles));
fprintf('Files updated: %d\n', successCount);
fprintf('Files with no changes: %d\n', noChangeCount);
fprintf('Errors: %d\n', errorCount);
fprintf('\nDone!\n');

%% Verification function - run this after the update
function verify_updates(labelDir)
    fprintf('\n=== VERIFICATION ===\n');
    matFiles = dir(fullfile(labelDir, "*.mat"));
    
    for k = 1:min(5, numel(matFiles))  % Check first 5 files
        matFile = fullfile(matFiles(k).folder, matFiles(k).name);
        
        try
            S = load(matFile, 'gTruth');
            gTruth = S.gTruth;
            
            fprintf('File %d: %s\n', k, matFiles(k).name);
            for dsIdx = 1:length(gTruth.DataSource)
                fprintf('  DataSource %d: "%s"\n', dsIdx, gTruth.DataSource(dsIdx).SourceName);
            end
            
        catch ME
            fprintf('File %d ERROR: %s\n', k, ME.message);
        end
    end
end