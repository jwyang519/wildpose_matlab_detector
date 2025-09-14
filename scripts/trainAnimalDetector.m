%% ==================== PointPillars Script with Simple Data Augmentation ====================

%% 0) Load multiple ground-truth files and merge
disp("🔄 Step 0: Loading all ground truth files...");

proj = currentProject;
root = proj.RootFolder;
matDir  = fullfile(root, "label/28-8-2025_30/");
pattern = "*.mat";

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

%% 1) Generate lidar training datastores
disp("🔄 Step 1: Generating training datastores...");
[ptds, blds] = lidarObjectDetectorTrainingData(gTruthAll);

%% 1.1) Apply MATLAB-based data augmentation
disp("🔄 Step 1.1: Setting up MATLAB data augmentation...");

pair = combine(ptds, blds);

% Optional: oversample objects (needs classNames and a sampled pool)
[ldsS,bdsS] = sampleLidarData(pair, classNames, 'MinPoints',20);
pair = transform(pair, @(x) pcBboxOversample(x, combine(ldsS,bdsS), classNames, [nCar nTruck]));

% Global aug in one place: the function must accept/return the combined row
pairAug = transform(pair, @(row) helperAugmentData(row));  % row is {ptCloud, boxes, labels}
cds = pairAug;      % hand this to trainPointPillarsObjectDetector

disp("✅ Step 1.1 complete: MATLAB augmentation applied successfully.");


%% 1.5) Extract full training table for anchor estimation
disp("🔄 Step 1.5: Extracting full training table for anchor/class info...");
nTruth = numel(gTruthAll);
samplingFactors = ones(1, nTruth);
fullData = lidarObjectDetectorTrainingData(gTruthAll, 'SamplingFactor', samplingFactors);
disp("✅ Step 1.5 complete: fullData rows: " + height(fullData))

classNames = fullData.Properties.VariableNames(2:end);

%% 2) Stratified split into train/val/test
disp("🔄 Step 2: Performing stratified splitting...");
classCols = classNames;
numClasses = numel(classCols);

trainIdx = false(height(fullData),1);
valIdx   = false(height(fullData),1);
testIdx  = false(height(fullData),1);

for i = 1:numClasses
    cls = classCols{i};
    idxCls = find(~cellfun(@isempty, fullData.(cls)));
    n = numel(idxCls);

    rng(1);
    idxShuf = idxCls(randperm(n));

    nTrain = floor(0.7 * n);
    nVal   = floor(0.15 * n);
    nTest  = n - (nTrain + nVal);

    trainIdx(idxShuf(1:nTrain)) = true;
    valIdx(idxShuf(nTrain+1 : nTrain+nVal)) = true;
    testIdx(idxShuf(nTrain+nVal+1 : end)) = true;
end

trainIdx = trainIdx & ~valIdx & ~testIdx;
valIdx   = valIdx & ~testIdx;

trainTable = fullData(trainIdx, :);
valTable   = fullData(valIdx, :);
testTable  = fullData(testIdx, :);

disp("✅ Step 2 complete: Stratified split sizes — Train: " + ...
    sum(trainIdx) + ", Val: " + sum(valIdx) + ", Test: " + sum(testIdx));

%% 3) Estimate pcRange from training point clouds
disp("🔄 Step 3: Estimating pcRange from training set...");
allRanges = zeros(height(trainTable),6);
for i = 1:height(trainTable)
    pc = pcread(trainTable.ptCloudFilename{i});
    xyz = reshape(pc.Location, [], 3);
    allRanges(i,:) = [ min(xyz(:,1)), max(xyz(:,1)), ...
                       min(xyz(:,2)), max(xyz(:,2)), ...
                       min(xyz(:,3)), max(xyz(:,3)) ];
end
pcRange = [min(allRanges(:,1)), max(allRanges(:,2)), ...
           min(allRanges(:,3)), max(allRanges(:,4)), ...
           min(allRanges(:,5)), max(allRanges(:,6))];
buffer = [0.5 0.5 0.5 0.5 0.2 0.2];
pcRange = pcRange + [-buffer(1:3), buffer(4:6)];
disp("✅ Step 3 complete: pcRange estimated.");

%% 4) Compute voxel size aligned to multiples of 8
xLen = pcRange(2) - pcRange(1);
yLen = pcRange(4) - pcRange(3);
xStep = 0.125; yStep = 0.125;
nX8 = ceil(xLen / xStep / 8) * 8;
nY8 = ceil(yLen / yStep / 8) * 8;
vx = xLen / nX8;
vy = yLen / nY8;
voxelSize = [vx vy];
disp("✅ Step 4 complete: Voxel size = [" + vx + ", " + vy + "]");

%% 5) Estimate anchor boxes
disp("🔄 Step 5: Estimating anchor boxes from train data...");
rawAnchors = calculateAnchorsPointPillars(trainTable(:,2:end));
anchorBoxes = rawAnchors(:);
for i = 1:numel(anchorBoxes)
    anchorsCi = anchorBoxes{i};
    if isempty(anchorsCi)
        warning("No anchors for class '%s'; inserting default placeholder anchors.", classNames{i});
        anchorsCi = [1 1 1 0 0; 1 1 1 0 pi/2];
    else
        uniqueYaws = unique(anchorsCi(:,5));
        if numel(uniqueYaws) < 2
            anchorsCi = [anchorsCi; anchorsCi(1,1:4), anchorsCi(1,5) + pi/2];
        end
    end
    anchorBoxes{i} = anchorsCi;
end
disp("✅ Step 5 complete: Anchor boxes ready.");

%% 6) Build untrained PointPillars detector
detector = pointPillarsObjectDetector(pcRange, classNames, anchorBoxes, ...
                                     VoxelSize=voxelSize);
disp("✅ Step 6 complete: Detector initialized.");

%% 7) Set training options with dynamic checkpoint path
disp("🔄 Step 7: Setting training options with dynamic checkpoint path...");

% Create timestamp for unique folder name
timestamp = string(datetime("now", "Format", "yyyy-MM-dd_HH-mm-ss"));
checkpointFolderName = sprintf('training_%s', timestamp);

% Create full checkpoint path
checkpointDir = fullfile(root, 'checkpoints', checkpointFolderName);

% Create the checkpoint directory
if ~exist(fullfile(root, 'checkpoints'), 'dir')
    mkdir(fullfile(root, 'checkpoints'));
    disp("Created parent checkpoints directory: " + fullfile(root, 'checkpoints'));
end

mkdir(checkpointDir);
disp("Created checkpoint directory: " + checkpointDir);

% Set training options
options = trainingOptions('adam', ...
    'Plots', "training-progress", ...
    'MaxEpochs', 50, ...
    'MiniBatchSize', 2, ...
    'ValidationFrequency', 100, ...
    'ValidationPatience', 8, ...
    'InitialLearnRate', 1e-4, ...
    'LearnRateSchedule', "piecewise", ...
    'LearnRateDropPeriod', 20, ...
    'LearnRateDropFactor', 0.8, ...
    'CheckpointPath', checkpointDir, ...
    'CheckpointFrequency', 25, ...
    'ResetInputNormalization', false, ...
    'Shuffle', "every-epoch");

disp("✅ Step 7 complete: Training options set.");
disp("📁 Checkpoint directory: " + checkpointDir);

%% 8) Train the model
disp("🚀 Step 8: Starting model training with augmentation...");
disp("📁 Training checkpoints will be saved to: " + checkpointDir);

[detector, info] = trainPointPillarsObjectDetector(cds_o, detector, options);
disp("✅ Step 8 complete: Model training finished.");

%% 9) Save model with timestamp
modelFilename = sprintf('trainedPointPillarsDetector_%s.mat', timestamp);
modelPath = fullfile(root, modelFilename);

save(modelPath, "detector", "info", "timestamp", "checkpointDir", "augParams");
disp("💾 Model saved to: " + modelPath);
disp("📁 Training checkpoints available at: " + checkpointDir);

%% ==================== AUGMENTATION FUNCTIONS ====================

function rowOut = helperAugmentData(rowIn)
% rowIn/rowOut: 1x3 cell {ptCloud, boxes, labels}
pc = rowIn{1}; boxes = rowIn{2}; labels = rowIn{3};

% rotate/scale/translate pc.Location and update boxes(:,1:3) and yaw
th = (rand*2-1) * (pi/4);    % example: ±45°
c=cos(th); s=sin(th); R=[c -s 0; s c 0; 0 0 1];

P = (R * pc.Location')';
if ~isempty(pc.Intensity), pc = pointCloud(P,'Intensity',pc.Intensity); else, pc = pointCloud(P); end

if ~isempty(boxes)
    XY = boxes(:,1:2);
    XY = (R(1:2,1:2)*XY')';     % rotate centers in XY
    boxes(:,1:2) = XY;
    if size(boxes,2) >= 7, boxes(:,7) = boxes(:,7) + th; end
end

rowOut = {pc, boxes, labels};  % same structure back
end


function anchors = calculateAnchorsPointPillars(lblTbl, useKmeans, k)
    % Create 3-D anchor boxes for PointPillars
    if nargin < 2, useKmeans = false; end
    if nargin < 3, k = 2; end

    nClasses = width(lblTbl);
    classNames = lblTbl.Properties.VariableNames;
    anchors = cell(1, nClasses);

    for ci = 1:nClasses
        colCell = lblTbl.(classNames{ci});
        allBoxes = cat(1, colCell{:});

        if isempty(allBoxes)
            warning("No boxes found for class '%s' – skipping.", classNames{ci});
            continue
        end

        dims = allBoxes(:, 4:6);
        zc = allBoxes(:, 3);

        if useKmeans
            [~, cent] = kmeans(dims, k, "Replicates", 5, "MaxIter", 500);
        else
            cent = mean(dims, 1);
        end

        if isvector(cent), cent = cent(:).'; end
        nA = size(cent, 1);
        meanZ = mean(zc);
        yaws = [0; pi/2];

        anchorsCi = zeros(nA * numel(yaws), 5);
        row = 1;
        for i = 1:nA
            for y = 1:numel(yaws)
                anchorsCi(row, :) = [cent(i, :) meanZ yaws(y)];
                row = row + 1;
            end
        end

        anchors{ci} = anchorsCi;
    end
end