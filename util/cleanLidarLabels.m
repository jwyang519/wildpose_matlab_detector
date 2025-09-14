function cleanLidarLabels(inputFolder, outputFolder, positionThreshold, maxKeep)
    % cleanLidarLabels - Remove repetitive bounding box labels from LiDAR gTruth files
    %
    % Inputs:
    %   inputFolder - Path to folder containing .mat files with gTruth objects
    %   outputFolder - Path to folder where cleaned files will be saved
    %   positionThreshold - Distance threshold for considering boxes as identical (default: 0.5)
    %   maxKeep - Maximum number of identical boxes to keep (default: 5)
    
    if nargin < 3
        positionThreshold = 0.5; % meters
    end
    if nargin < 4
        maxKeep = 5;
    end
    
    % Create output folder if it doesn't exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Get all .mat files in input folder
    matFiles = dir(fullfile(inputFolder, '*.mat'));
    
    fprintf('Found %d .mat files to process...\n', length(matFiles));
    
    for fileIdx = 1:length(matFiles)
        filename = matFiles(fileIdx).name;
        filepath = fullfile(inputFolder, filename);
        
        fprintf('Processing %s...\n', filename);
        
        % Load the gTruth object
        data = load(filepath);
        fieldNames = fieldnames(data);
        
        % Find the gTruth object (try multiple possible class names and field names)
        gTruthObj = [];
        gTruthFieldName = '';
        
        % First, try common field names
        commonFieldNames = {'gTruth', 'gt', 'groundTruth', 'labels'};
        for i = 1:length(commonFieldNames)
            if isfield(data, commonFieldNames{i})
                candidate = data.(commonFieldNames{i});
                if isgTruthObject(candidate)
                    gTruthObj = candidate;
                    gTruthFieldName = commonFieldNames{i};
                    break;
                end
            end
        end
        
        % If not found, check all fields
        if isempty(gTruthObj)
            for i = 1:length(fieldNames)
                candidate = data.(fieldNames{i});
                if isgTruthObject(candidate)
                    gTruthObj = candidate;
                    gTruthFieldName = fieldNames{i};
                    break;
                end
            end
        end
        
        if isempty(gTruthObj)
            fprintf('Warning: No gTruth object found in %s\n', filename);
            fprintf('Available fields: %s\n', strjoin(fieldNames, ', '));
            % Try to show the class of each field for debugging
            for i = 1:length(fieldNames)
                fprintf('  %s: %s\n', fieldNames{i}, class(data.(fieldNames{i})));
            end
            continue;
        end
        
        % Clean the labels
        cleanedGTruth = cleanRepetitiveBoxes(gTruthObj, positionThreshold, maxKeep);
        
        % Save cleaned version
        outputPath = fullfile(outputFolder, filename);
        data.(gTruthFieldName) = cleanedGTruth;
        save(outputPath, '-struct', 'data');
        
        fprintf('Cleaned %s saved to %s\n', filename, outputPath);
    end
    
    fprintf('All files processed!\n');
end

function cleanedGTruth = cleanRepetitiveBoxes(gTruth, positionThreshold, maxKeep)
    % Extract all bounding boxes from all frames
    allBoxes = [];
    frameIndices = [];
    
    labelNames = gTruth.LabelDefinitions.Name;
    
    % Collect all bounding boxes across all frames
    for frameIdx = 1:height(gTruth.LabelData)
        for labelIdx = 1:length(labelNames)
            labelName = labelNames{labelIdx};
            
            if ismember(labelName, gTruth.LabelData.Properties.VariableNames)
                boxes = gTruth.LabelData.(labelName){frameIdx};
                
                if ~isempty(boxes)
                    % Store box info: [x, y, z, length, width, height, yaw, class_id, frame, original_row_in_frame]
                    numBoxes = size(boxes, 1);
                    classId = labelIdx * ones(numBoxes, 1);
                    frameId = frameIdx * ones(numBoxes, 1);
                    rowIds = (1:numBoxes)';
                    
                    boxInfo = [boxes, classId, frameId, rowIds];
                    allBoxes = [allBoxes; boxInfo];
                    frameIndices = [frameIndices; repmat(frameIdx, numBoxes, 1)];
                end
            end
        end
    end
    
    if isempty(allBoxes)
        fprintf('No bounding boxes found in this file.\n');
        cleanedGTruth = gTruth;
        return;
    end
    
    fprintf('Total boxes before cleaning: %d\n', size(allBoxes, 1));
    
    % Group similar boxes
    boxGroups = groupSimilarBoxes(allBoxes, positionThreshold);
    
    % Determine which boxes to keep
    boxesToRemove = [];
    totalRemoved = 0;
    
    for groupIdx = 1:length(boxGroups)
        group = boxGroups{groupIdx};
        
        if length(group) > maxKeep
            % Sort by frame index to keep boxes spread across time
            [~, sortIdx] = sort(allBoxes(group, end-1)); % Sort by frame index
            group = group(sortIdx);
            
            % Keep first maxKeep boxes, mark rest for removal
            boxesToRemove = [boxesToRemove; group(maxKeep+1:end)];
            totalRemoved = totalRemoved + length(group) - maxKeep;
        end
    end
    
    fprintf('Removing %d repetitive boxes...\n', totalRemoved);
    
    % Create new label data with cleaned boxes
    newLabelData = gTruth.LabelData;
    
    % Remove boxes from each frame
    for frameIdx = 1:height(newLabelData)
        for labelIdx = 1:length(labelNames)
            labelName = labelNames{labelIdx};
            
            if ismember(labelName, newLabelData.Properties.VariableNames)
                boxes = newLabelData.(labelName){frameIdx};
                
                if ~isempty(boxes)
                    % Find boxes to remove from this frame/label combination
                    frameBoxIndices = find(allBoxes(:, end-1) == frameIdx & allBoxes(:, end-2) == labelIdx);
                    boxesToRemoveThisFrame = [];
                    
                    for i = 1:length(frameBoxIndices)
                        globalBoxIdx = frameBoxIndices(i);
                        if ismember(globalBoxIdx, boxesToRemove)
                            localBoxIdx = allBoxes(globalBoxIdx, end); % original row in frame
                            boxesToRemoveThisFrame = [boxesToRemoveThisFrame; localBoxIdx];
                        end
                    end
                    
                    % Remove the boxes
                    if ~isempty(boxesToRemoveThisFrame)
                        boxes(boxesToRemoveThisFrame, :) = [];
                        newLabelData.(labelName){frameIdx} = boxes;
                    end
                end
            end
        end
    end
    
    % Create new gTruth object with cleaned data
    try
        % Try to create a new groundTruthLidar object
        cleanedGTruth = groundTruthLidar(gTruth.DataSource, gTruth.LabelDefinitions, newLabelData);
    catch ME
        fprintf('Error creating new groundTruthLidar object: %s\n', ME.message);
        % If that fails, try alternative approaches
        try
            % Method 2: Try copying and then replacing via constructor parameters
            cleanedGTruth = groundTruthLidar(gTruth.DataSource, gTruth.LabelDefinitions, newLabelData);
        catch ME2
            fprintf('Alternative method also failed: %s\n', ME2.message);
            % Method 3: Create a basic groundTruth object if groundTruthLidar fails
            try
                cleanedGTruth = groundTruth(gTruth.DataSource, gTruth.LabelDefinitions, newLabelData);
            catch ME3
                fprintf('All methods failed. Returning original gTruth: %s\n', ME3.message);
                cleanedGTruth = gTruth;
            end
        end
    end
    
    fprintf('Cleaning complete!\n');
end

function boxGroups = groupSimilarBoxes(allBoxes, threshold)
    % Group boxes that are similar in position and class
    % allBoxes format: [x, y, z, length, width, height, yaw, class_id, frame, row_id]
    
    numBoxes = size(allBoxes, 1);
    visited = false(numBoxes, 1);
    boxGroups = {};
    
    for i = 1:numBoxes
        if visited(i)
            continue;
        end
        
        currentGroup = i;
        visited(i) = true;
        
        % Find all similar boxes
        for j = i+1:numBoxes
            if visited(j)
                continue;
            end
            
            % Check if boxes are similar
            if areBoxesSimilar(allBoxes(i, :), allBoxes(j, :), threshold)
                currentGroup = [currentGroup; j];
                visited(j) = true;
            end
        end
        
        % Only consider groups with more than 1 box
        if length(currentGroup) > 1
            boxGroups{end+1} = currentGroup;
        end
    end
end

function similar = areBoxesSimilar(box1, box2, threshold)
    % Check if two boxes are similar based on class and position
    
    % Must be same class
    if box1(end-2) ~= box2(end-2)  % class_id comparison
        similar = false;
        return;
    end
    
    % Calculate Euclidean distance between centers
    pos1 = box1(1:3);  % [x, y, z]
    pos2 = box2(1:3);
    distance = norm(pos1 - pos2);
    
    % Check size similarity (optional - you can adjust or remove this)
    size1 = box1(4:6);  % [length, width, height]
    size2 = box2(4:6);
    sizeRatio = max(size1./size2, size2./size1);
    sizeSimilar = all(sizeRatio < 1.5);  % Allow 50% size difference
    
    similar = distance < threshold && sizeSimilar;
end

function isGTruth = isgTruthObject(obj)
    % Check if an object is a gTruth object using multiple criteria
    
    isGTruth = false;
    
    % Check class names (try different possible class names)
    possibleClasses = {'groundTruth', 'gTruth', 'GroundTruth'};
    objClass = class(obj);
    
    if any(strcmp(objClass, possibleClasses))
        isGTruth = true;
        return;
    end
    
    % Check if it has gTruth-like properties
    if isstruct(obj) || isobject(obj)
        try
            % Check for typical gTruth properties
            hasLabelData = isprop(obj, 'LabelData') || isfield(obj, 'LabelData');
            hasLabelDefinitions = isprop(obj, 'LabelDefinitions') || isfield(obj, 'LabelDefinitions');
            
            if hasLabelData && hasLabelDefinitions
                isGTruth = true;
                return;
            end
        catch
            % If property checking fails, continue
        end
    end
end

% Example usage:
% cleanLidarLabels('input_labels/', 'cleaned_labels/', 0.5, 5);