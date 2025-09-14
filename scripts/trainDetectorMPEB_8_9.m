%% ==================== Full PointPillars Script with Datastore ====================

%% 0) Load multiple ground-truth files and merge
disp("🔄 Step 0: Loading all ground truth files...");
proj = matlab.project.rootProject;
rootDir = proj.RootFolder;
matDir  = fullfile(rootDir, 'labels','8-9-2025_downsampled');
pattern = "*.mat";
outputFolder =fullfile(rootDir, 'TrainingOutput_08-09-2025');

files = dir(fullfile(matDir, pattern));
% Filter out macOS sidecars / hidden files
mask  = ~startsWith({files.name}, '.') & ~startsWith({files.name}, '._') ...
        & ~strcmp({files.name}, '.DS_Store');
files = files(mask);

gTruthAll = groundTruthLidar.empty;

for k = 1:numel(files)
    inFile = fullfile(files(k).folder, files(k).name);
    S = load(inFile);

    % Find the groundTruthLidar variable inside the MAT
    g = [];
    fns = fieldnames(S);
    for j = 1:numel(fns)
        if isa(S.(fns{j}), 'groundTruthLidar')
            g = S.(fns{j});
            break;
        end
    end
    if isempty(g)
        warning('[SKIP] %s (no groundTruthLidar found)', files(k).name);
        continue;
    end

    gTruthAll(end+1,1) = g; %#ok<SAGROW>
end

disp("✅ Step 0 complete: " + numel(gTruthAll) + " gTruth files loaded.");
%% 1) Manually setting pcRange and compute voxel size aligned to multiples of 8
disp("🔄 Step 1:Setting pcRange and voxel sizes...");
xMin = 0.0;     % Minimum value along X-axis.
yMin = -20;  % Minimum value along Y-axis.
zMin = -5.0;    % Minimum value along Z-axis.
xMax = 200;   % Maximum value along X-axis.
yMax = 20;   % Maximum value along Y-axis.
zMax = 15.0;     % Maximum value along Z-axis.

% Define point cloud parameters.
pcRange = [xMin xMax yMin yMax zMin zMax];

xLen = pcRange(2) - pcRange(1);
yLen = pcRange(4) - pcRange(3);
xStep = 0.15; yStep = 0.15;
nX8 = ceil(xLen / xStep / 8) * 8;
nY8 = ceil(yLen / yStep / 8) * 8;
vx = xLen / nX8;
vy = yLen / nY8;
voxelSize = [vx vy];
disp("✅ Step 1 complete: Voxel size = [" + vx + ", " + vy + "]");

%% 2) Generate lidar training datastores
disp("🔄 Step 2: Generating training datastores...");
nTruth = numel(gTruthAll);
samplingFactors = ones(1, nTruth);
[pcds, ~] = lidarObjectDetectorTrainingData(gTruthAll, 'SamplingFactor', samplingFactors);
fullData = lidarObjectDetectorTrainingData(gTruthAll, 'SamplingFactor', samplingFactors);
n = height(fullData);
disp("✅ Step 2 complete: fullData rows: " + n)

classNames = fullData.Properties.VariableNames(2:end);

disp("🔄 Cropping PCDs front view...");
boxLabels = fullData(:,2:5);
[croppedPointCloudObj,processedLabels] = cropFrontViewFromLidarData(...
    pcds,boxLabels,pcRange);
disp("✅ PCD objects cropped.");

%% 3) Stratified split into train/val/test from fullData table
disp("🔄 Step 3: Spliting data 70/30...");
rng(1);
shuffledIndices = randperm(size(processedLabels,1));
idx = floor(0.7 * length(shuffledIndices));

trainData = croppedPointCloudObj(shuffledIndices(1:idx),:);
testData = croppedPointCloudObj(shuffledIndices(idx+1:end),:);

trainLabels = processedLabels(shuffledIndices(1:idx),:);
testLabels = processedLabels(shuffledIndices(idx+1:end),:);

writeFiles = true;
dataLocation = fullfile(outputFolder,'InputData');
[trainData,trainLabels] = saveptCldToPCD(trainData,trainLabels,...
    dataLocation,writeFiles);

lds = fileDatastore(dataLocation,'ReadFcn',@(x) pcread(x));
bds = boxLabelDatastore(trainLabels);
cds = combine(lds,bds);
disp("✅ Step 3: Spliting done, file datastore saved.");

%% 4) Data augmentation
disp("Skipping data augmentation for now...");

%% 5) Estimate anchor boxes
disp("🔄 Step 5: Estimating anchor boxes from train data...");
anchorBoxes = calculateAnchorsPointPillars(trainLabels);
% rawAnchors = calculateAnchorsPointPillars(trainTable(:,2:end));
% anchorBoxes = rawAnchors(:);
% for i = 1:numel(anchorBoxes)
%     anchorsCi = anchorBoxes{i};
%     if isempty(anchorsCi)
%         warning("No anchors for class '%s'; inserting default placeholder anchors.", classNames{i});
%         anchorsCi = [1 1 1 0 0; 1 1 1 0 pi/2];
%     else
%         uniqueYaws = unique(anchorsCi(:,5));
%         if numel(uniqueYaws) < 2
%             anchorsCi = [anchorsCi; anchorsCi(1,1:4), anchorsCi(1,5) + pi/2];
%         end
%     end
%     anchorBoxes{i} = anchorsCi;
% end
disp("✅ Step 5 complete: Anchor boxes ready.");

%% 6) Build untrained PointPillars detector
detector = pointPillarsObjectDetector(pcRange, classNames, anchorBoxes, ...
                                     VoxelSize=voxelSize);
disp("✅ Step 6 complete: Detector initialized.");

%% 7) Set training options
% Create timestamp for unique folder name
timestamp = string(datetime("now", "Format", "yyyy-MM-dd_HH-mm-ss"));
checkpointFolderName = sprintf('training_%s', timestamp);

% Create full checkpoint path
checkpointDir = fullfile(outputFolder, 'checkpoints', checkpointFolderName);

% Create the checkpoint directory
if ~exist(fullfile(outputFolder, 'checkpoints'), 'dir')
    mkdir(fullfile(outputFolder, 'checkpoints'));
    disp("Created parent checkpoints directory: " + fullfile(outputFolder, 'checkpoints'));
end

mkdir(checkpointDir);
disp("Created checkpoint directory: " + checkpointDir);
options = trainingOptions('adam', ...
    Plots='training-progress', ...
    MaxEpochs=100, ...
    MiniBatchSize=1, ...
    GradientDecayFactor = 0.9,...
    SquaredGradientDecayFactor = 0.999,...
    InitialLearnRate= 2e-4 , ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropPeriod=15, ...
    LearnRateDropFactor=0.8, ...
    CheckpointPath=checkpointDir, ...
    CheckpointFrequency=5, ...
    CheckpointFrequencyUnit = 'epoch', ...
    ResetInputNormalization=false, ...
    executionEnvironment="auto");
disp("✅ Step 7 complete: Training options set.");

%% 8) Train using the combined datastore
disp("🚀 Step 8: Starting model training...");
[detector, info] = trainPointPillarsObjectDetector(cds, detector, options);
disp("✅ Step 8 complete: Model training finished.");

%% 9) Save model
save("trainedPointPillarsDetector.mat", "detector", "info");
disp("💾 Model saved to trainedPointPillarsDetector.mat");

function [croppedPointCloudObj, processedLabels] = cropFrontViewFromLidarData(lidarData, boxLabels, pcRange)
% This function crops the front view from the input full-view point cloud
% and also processes the corresponding box labels according to the 
% specified grid parameters.

    xmin = pcRange(1,1);
    xmax = pcRange(1,2);
    ymin = pcRange(1,3);
    ymax = pcRange(1,4);
    zmin = pcRange(1,5);
    zmax = pcRange(1,6);
    
    tmpStr = '';
    numFiles = size(boxLabels,1);
    
    processedLabels = cell(size(boxLabels));
    croppedPointCloudObj = cell(size(numFiles));

    % Get the classes from the ground truth labels.
    classNames = boxLabels.Properties.VariableNames;
    
    for i = 1:numFiles

        ptCloud = read(lidarData);            
        groundTruth = boxLabels(i,:);
        
        % Set the limits for the point cloud.
        if(numel(size(ptCloud.Location)) == 3)
            % Organized point cloud
            [x,y] = find( ptCloud.Location(:,:,1) < xmax ...
                                & ptCloud.Location(:,:,1) > xmin ...
                                & ptCloud.Location(:,:,2) < ymax ...
                                & ptCloud.Location(:,:,2) > ymin ...
                                & ptCloud.Location(:,:,3) < zmax ...
                                & ptCloud.Location(:,:,3) > zmin);    
            ptCloud = select(ptCloud, x, y, 'OutputSize', 'full'); 
        else
            % Unorganized point cloud
            pos = find( ptCloud.Location(:,1) < xmax ...
                & ptCloud.Location(:,1) > xmin ...
                & ptCloud.Location(:,2) < ymax ...
                & ptCloud.Location(:,2) > ymin ...
                & ptCloud.Location(:,3) < zmax ...
                & ptCloud.Location(:,3) > zmin);    
            ptCloud = select(ptCloud, pos, 'OutputSize', 'full'); 
        end
        
        processedData = removeInvalidPoints(ptCloud);
         
        for ii = 1:numel(classNames)

            labels = groundTruth(1,classNames{ii}).Variables;
            if(iscell(labels))
                labels = labels{1};
            end
            if ~isempty(labels)

                % Get the label indices that are in the selected RoI.
                labelsIndices = labels(:,1) > xmin ...
                            & labels(:,1) < xmax ...
                            & labels(:,2) > ymin...
                            & labels(:,2) < ymax ...
                            & labels(:,4) > 0 ...
                            & labels(:,5) > 0 ...
                            & labels(:,6) > 0;
                labels = labels(labelsIndices,:);

                if ~isempty(labels)
                    % Find the number of points inside each ground truth
                    % label.
                    numPoints = arrayfun(@(x)(findPointsInModel(cuboidModel(labels(x,:)),processedData)),...
                                (1:size(labels,1)).','UniformOutput',false);

                    posLabels = cellfun(@(x)(length(x) > 50), numPoints);
                    labels = labels(posLabels,:);
                end
            end
            processedLabels{i,ii} = labels;
        end
        croppedPointCloudObj{i,1} = processedData;
    end
    
    % Print completion message when done.
    msg = sprintf('Processing data 100%% complete');
    fprintf(1,'%s',[tmpStr, msg]);

    processedLabels = cell2table(processedLabels);
    numClasses = size(processedLabels,2);
    for j = 1:numClasses
        processedLabels.Properties.VariableNames{j} = classNames{j};
    end

end

function [ptCld, ptLabels] = saveptCldToPCD(ptCld, ptLabels, dataLocation, writeFiles)
% This function saves the required point clouds in the specified location 

    if ~exist(dataLocation, 'dir')
        mkdir(dataLocation)
    end
    
    tmpStr = '';
    numFiles = size(ptLabels,1);
    ind = [];
    
    for i = 1:numFiles
        processedData = ptCld{i,1};
        
        % Skip if the processed point cloud is empty
        if(isempty(processedData.Location))
            ind = [ind, i];
            continue;
        end
        
        if(writeFiles)
            pcFilePath = fullfile(dataLocation, sprintf('%06d.pcd',i));
            pcwrite(processedData, pcFilePath);
        end
      
        % Display progress after 300 files on screen.
        if ~mod(i,300)
            msg = sprintf('Processing data %3.2f%% complete', (i/numFiles)*100.0);
            fprintf(1,'%s',[tmpStr, msg]);
            tmpStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
    
    ptCld(ind,:) = [];
    ptLabels(ind,:) = [];
    
    % Print completion message when done.
    msg = sprintf('Processing data 100%% complete');
    fprintf(1,'%s',[tmpStr, msg]);
end

function anchors = calculateAnchorsPointPillars(labels)
% This function calculates the anchor boxes from the box labels in the
% training set.

    anchors = [];
    classNames = labels.Properties.VariableNames;
    
    % Calculate the anchors for each class label.
    for ii = 1:numel(classNames)
        bboxCells = table2array(labels(:,ii));
        lwhValues = [];
        
        % Accumulate the lengths, widths, heights from the ground truth
        % labels.
        for i = 1 : height(bboxCells)
            if(~isempty(bboxCells{i}))
                lwhValues = [lwhValues; bboxCells{i}(:, 4:6)];
            end
        end
        
        % Calculate the mean for each. 
        meanVal = mean(lwhValues, 1);
        
        % With the obtained mean values, create two anchors with two 
        % yaw angles, 0 and 90.
%         classAnchors = [{num2cell([meanVal, -1.78, 0])}, {num2cell([meanVal, -1.78, pi/2])}];
        classAnchors = [[meanVal, -1.78, 0]; [meanVal, -1.78, pi/2]];
        
        anchors = [anchors; {classAnchors}];
    end
end