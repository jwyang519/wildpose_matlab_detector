%% ==================== Full PointPillars Script with Datastore ====================

%% 0) Load multiple ground-truth files and merge
disp("🔄 Step 0: Loading all ground truth files...");

rootDir = "/Volumes/T7_Shield/WildPose/WildPose_Proj";
matDir  = "/Volumes/T7_Shield/WildPose/WildPose_Proj/label/28-8-2025_30";
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
cds = combine(ptds, blds);
disp("✅ Step 1 complete: Combined datastore created.");
% the above data was working when running in the HPC, just needs the
% following sanitisation

% --- ADD THIS right after you create `cds` (and similarly for `valPair`) ---
epsilon = 1e-3;   % minimal positive size for [w,l,h]

sanitizeFcn = @(row) iSanitizeRow(row, epsilon);

% For training:
cds = transform(cds, sanitizeFcn);

% (Optional) quick pass to ensure no bad rows remain
try
    reset(cds);
    while hasdata(cds)
        read(cds);
    end
    reset(cds);
catch ME
    error("Sanitized datastore still contains a bad row: %s", ME.message);
end
%% 1.5) Also generate fullData table for analysis/anchor estimation
disp("🔄 Step 1.5: Extracting full training table for anchor/class info...");
nTruth = numel(gTruthAll);
samplingFactors = ones(1, nTruth);
fullData = lidarObjectDetectorTrainingData(gTruthAll, 'SamplingFactor', samplingFactors);
disp("✅ Step 1.5 complete: fullData rows: " + height(fullData))

classNames = fullData.Properties.VariableNames(2:end);

%% 2) Stratified split into train/val/test from fullData table
disp("🔄 Step 2: Performing stratified splitting...");

% 1) Create full dataset to get total number of rows
fullData = lidarObjectDetectorTrainingData(gTruthAll, 'SamplingFactor', ones(1,numel(gTruthAll)));
n = height(fullData);

% 2) Row-level split using 3:1:1 ratio (period of 5)
period = 5; 
r = mod(0:n-1, period);
trainMask = r < 3;      % rows 0,1,2 -> 60% training
valMask = r == 3;       % row 3 -> 20% validation  
testMask = r == 4;      % row 4 -> 20% test

% 3) Build row→file mapping (same as original)
counts = zeros(numel(gTruthAll),1);
for gi = 1:numel(gTruthAll)
    Ti = lidarObjectDetectorTrainingData(gTruthAll(gi), 'SamplingFactor', 1);
    counts(gi) = height(Ti);
end
row2gt = repelem((1:numel(gTruthAll))', counts);

% 4) For each split, collect the corresponding ground truth files
% But keep ALL files that have ANY rows in each split (row-level, not file-level)
gtIdxTrain = unique(row2gt(trainMask));
gtIdxVal = unique(row2gt(valMask));  
gtIdxTest = unique(row2gt(testMask));

% 5) Create ground truth objects for each split
gTruthTrain = gTruthAll(gtIdxTrain);
gTruthVal = gTruthAll(gtIdxVal);
gTruthTest = gTruthAll(gtIdxTest);

% 6) Display results
disp("✅ Step 1 complete: Row-level split sizes — Train: " + ...
     numel(gTruthTrain) + ", Val: " + numel(gTruthVal) + ", Test: " + numel(gTruthTest));

% 7) Show actual row counts for verification
trainRows = sum(trainMask);
valRows = sum(valMask);
testRows = sum(testMask);

disp("📊 Actual row distribution — Train: " + trainRows + ...
     " (" + sprintf('%.1f%%', trainRows/n*100) + "), Val: " + valRows + ...
     " (" + sprintf('%.1f%%', valRows/n*100) + "), Test: " + testRows + ...
     " (" + sprintf('%.1f%%', testRows/n*100) + ")");

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
xStep = 0.25; yStep = 0.25;
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

%% 7) Set training options
% Create timestamp for unique folder name
timestamp = string(datetime("now", "Format", "yyyy-MM-dd_HH-mm-ss"));
checkpointFolderName = sprintf('training_%s', timestamp);

% Create full checkpoint path
checkpointDir = fullfile(rootDir, 'checkpoints', checkpointFolderName);

% Create the checkpoint directory
if ~exist(fullfile(rootDir, 'checkpoints'), 'dir')
    mkdir(fullfile(rootDir, 'checkpoints'));
    disp("Created parent checkpoints directory: " + fullfile(rootDir, 'checkpoints'));
end

mkdir(checkpointDir);
disp("Created checkpoint directory: " + checkpointDir);
options = trainingOptions('adam', ...
    Plots='none', ...
    MaxEpochs=100, ...
    MiniBatchSize=1, ...
    InitialLearnRate= 1e-4 , ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropPeriod=15, ...
    LearnRateDropFactor=0.8, ...
    CheckpointPath=checkpointDir, ...
    CheckpointFrequency=5, ...
    ResetInputNormalization=false, ...
    BatchNormalizationStatistics="moving");
disp("✅ Step 7 complete: Training options set.");

%% 8) Train using the combined datastore
disp("🚀 Step 8: Starting model training...");
[detector, info] = trainPointPillarsObjectDetector(cds, detector, options);
disp("✅ Step 8 complete: Model training finished.");

%% 9) Save model
save("trainedPointPillarsDetector.mat", "detector", "info");
disp("💾 Model saved to trainedPointPillarsDetector.mat");


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


function rowOut = iSanitizeRow(rowIn, epsilon)
% rowIn/rowOut: {ptCloud, boxes, labels}
pc    = rowIn{1};
boxes = rowIn{2};
labs  = rowIn{3};

% Empty or no boxes? pass through.
if isempty(boxes)
    rowOut = rowIn; 
    return;
end

% 3D box format is [x y z w l h yaw] (+ optional 8/9 cols). 
% Clip dims to be strictly positive; also clean NaN/Inf.
dims = boxes(:,4:6);
badNum = any(~isfinite(boxes),2);
dims = max(dims, epsilon);
boxes(:,4:6) = dims;

% Optionally normalize yaw (keeps it tidy; not strictly required)
if size(boxes,2) >= 7 && all(isfinite(boxes(:,7)))
    boxes(:,7) = wrapToPi(boxes(:,7));
end

% Drop boxes that are still invalid (e.g., NaNs in centers)
stillBad = badNum | any(~isfinite(boxes(:,1:7)),2) | any(boxes(:,4:6) <= 0,2);
if any(stillBad)
    boxes(stillBad,:) = [];
    labs(stillBad,:)  = [];
end

% If no boxes survive, return the cloud with empty labels (some trainers accept; PointPillars does)
rowOut = {pc, boxes, labs};
end