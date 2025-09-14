%% ==================== PointPillars Script with Simple Data Augmentation ====================

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

% Derive class names early so augmentation/oversampling can use them
tmpTbl = lidarObjectDetectorTrainingData(gTruthAll, 'SamplingFactor', ones(1, numel(gTruthAll)));
classNames = tmpTbl.Properties.VariableNames(2:end);

disp("✅ Step 0 complete: " + numel(gTruthAll) + " gTruth files loaded.");

% 0.5) Memory Check
try
    g = gpuDevice;
    fprintf('GPU: %s | Total %.1f GB | Free %.1f GB\n', g.Name, g.TotalMemory/1e9, g.AvailableMemory/1e9);
catch ME
    fprintf('gpuDevice() failed: %s\n', ME.message);
end

%% 1) Simple stratified split into train/val/test
disp("🔄 Step 1: Performing row-level periodic splitting...");

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
%% 1.a) DEBUG: class presence per split (ROW-LEVEL, using fullData masks)
fprintf('\n==== CLASS PRESENCE PER SPLIT (ROW-LEVEL via masks) ====\n');

% Build per-split row tables directly from fullData using the masks
Ttrain = fullData(trainMask, :);
Tval   = fullData(valMask,   :);
Ttest  = fullData(testMask,  :);

% Union of class names across all three (skip the first column)
classesTrain = string(Ttrain.Properties.VariableNames(2:end));
classesVal   = string(Tval.Properties.VariableNames(2:end));
classesTest  = string(Ttest.Properties.VariableNames(2:end));
allClasses   = unique([classesTrain, classesVal, classesTest]);

% Helper: count ROWS (frames) that contain at least one instance of a class
countRows = @(T, cls) arrayfun(@(k) ...
    (ismember(cls(k), string(T.Properties.VariableNames))) * ...
    sum(~cellfun(@isempty, T.(cls(k)))) , 1:numel(cls));

cnt.train = countRows(Ttrain, allClasses);
cnt.val   = countRows(Tval,   allClasses);
cnt.test  = countRows(Ttest,  allClasses);

% Print counts + totals
fprintf('%-15s %8s %8s %8s %8s\n', 'Class', 'Train', 'Val', 'Test', 'Total');
for i = 1:numel(allClasses)
    tot = cnt.train(i) + cnt.val(i) + cnt.test(i);
    fprintf('%-15s %8d %8d %8d %8d\n', allClasses(i), cnt.train(i), cnt.val(i), cnt.test(i), tot);
end

% Print percentages within each class
fprintf('\n%-15s %8s %8s %8s\n', 'Class', 'Train%', 'Val%', 'Test%');
for i = 1:numel(allClasses)
    tot = cnt.train(i) + cnt.val(i) + cnt.test(i);
    if tot > 0
        fprintf('%-15s %7.1f%% %7.1f%% %7.1f%%\n', allClasses(i), ...
            100*cnt.train(i)/tot, 100*cnt.val(i)/tot, 100*cnt.test(i)/tot);
    else
        fprintf('%-15s %8s %8s %8s\n', allClasses(i), 'N/A', 'N/A', 'N/A');
    end
end

% Warnings
missingInTrain = allClasses(cnt.train==0 & (cnt.val>0 | cnt.test>0));
if ~isempty(missingInTrain)
    warning('Classes present in VAL/TEST but missing in TRAIN: %s', strjoin(cellstr(missingInTrain), ', '));
end

absentAll = allClasses(cnt.train==0 & cnt.val==0 & cnt.test==0);
if ~isempty(absentAll)
    warning('Classes absent in ALL splits: %s', strjoin(cellstr(absentAll), ', '));
end

fprintf('===============================================\n\n');
%% 2) Generate training datastores from train split
disp("🔄 Step 2: Generating training datastores...");
[ptds, blds] = lidarObjectDetectorTrainingData(gTruthTrain);

%% 2.1) Apply MATLAB-based data augmentation to training data
disp("🔄 Step 2.1: Setting up MATLAB data augmentation...");

% Start from raw paired datastore
pair = combine(ptds, blds);

% ----- Test basic datastore first -----
try
    disp("Testing basic datastore read...");
    testRead = read(pair);
    reset(pair);
    disp("✅ Basic datastore read successful");
    disp("Data structure: " + class(testRead) + ", Size: " + numel(testRead));
    if iscell(testRead)
        for i = 1:numel(testRead)
            disp("  Cell " + i + ": " + class(testRead{i}));
        end
    end
catch ME
    error("Basic datastore read failed: %s", ME.message);
end

% ----- Skip oversampling for now - focus on geometric augmentation -----
disp("Skipping oversampling - using geometric augmentation only");

% ----- Test simple augmentation first -----
try
    disp("Testing geometric augmentation function...");
    testRow = read(pair);
    reset(pair);
    augmentedRow = helperAugmentData(testRow);
    disp("✅ Geometric augmentation function test successful");
catch ME
    error("Geometric augmentation function test failed: %s", ME.message);
end

% ----- Apply geometric augmentation only -----
cds = transform(pair, @(row) sanitizeRow(helperAugmentData(row)));

% (Optional) quick sanity check
try
    disp("Testing final augmented datastore...");
    t = read(cds); %#ok<NASGU>
    reset(cds);
    disp("✅ Final augmented datastore read successful");
catch ME
    error("Augmented datastore failed basic read: %s", ME.message);
end

disp("✅ Step 2.1 complete: MATLAB augmentation applied successfully.");

%% 2.2) Create validation datastores (NO augmentation)
disp("🔄 Step 2.2: Creating validation datastores...");
[ptdsVal, bldsVal] = lidarObjectDetectorTrainingData(gTruthVal);
valPairRaw = combine(ptdsVal, bldsVal);
valPair = transform(valPairRaw, @sanitizeRow);

% Optional quick test
vt = read(valPair); reset(valPair);
disp("✅ Step 2.2 complete: Validation datastores created (no augmentation).");

%% 3) Extract training table for anchor estimation and pcRange
disp("🔄 Step 3: Extracting training table for anchor/pcRange estimation...");
trainTable = lidarObjectDetectorTrainingData(gTruthTrain, 'SamplingFactor', ones(1, numel(gTruthTrain)));
disp("✅ Step 3 complete: trainTable rows: " + height(trainTable));

%% 4) Estimate pcRange from training point clouds
disp("🔄 Step 4: Estimating pcRange from training set...");
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
disp("✅ Step 4 complete: pcRange estimated.");

%% 4.1) Simple forward cap (keep Y/Z full span)
% disp("🔄 Step 4.1: Capping pcRange forward distance...");
% 
% oldRange = pcRange;
% 
% % Set forward cap 
% xFwdCap = 250;
% 
% % Keep original Y/Z ranges (to avoid cutting side/tall objects)
% xLo = oldRange(1);
% xHi = min(oldRange(2), xFwdCap);
% yLo = oldRange(3);
% yHi = oldRange(4);
% zLo = oldRange(5);
% zHi = oldRange(6);
% 
% pcRange = [xLo xHi yLo yHi zLo zHi];
% 
% fprintf("pcRange capped from [%.1f %.1f | %.1f %.1f | %.1f %.1f] -> [%.1f %.1f | %.1f %.1f | %.1f %.1f]\n", ...
%     oldRange(1),oldRange(2),oldRange(3),oldRange(4),oldRange(5),oldRange(6), ...
%     pcRange(1),pcRange(2),pcRange(3),pcRange(4),pcRange(5),pcRange(6));
% 
% % --- Coverage check: warn if any class lost all its boxes ---
% classes = trainTable.Properties.VariableNames(2:end);
% lost = [];
% for ci = 1:numel(classes)
%     col = trainTable.(classes{ci});
%     keepCount = 0;
%     for r = 1:numel(col)
%         B = col{r};
%         if isempty(B), continue; end
%         xc = B(:,1);  % box centers (x)
%         in = xc >= pcRange(1) & xc <= pcRange(2);
%         keepCount = keepCount + nnz(in);
%     end
%     if keepCount == 0
%         lost(end+1) = ci; %#ok<AGROW>
%     end
% end
% 
% if ~isempty(lost)
%     fprintf("⚠️ WARNING: These classes lost all boxes in capped range: %s\n", ...
%         strjoin(classes(lost), ', '));
% else
%     disp("✅ All classes still have labeled boxes within capped range.");
% end
%% 5) Compute voxel size aligned to multiples of 8
xLen = pcRange(2) - pcRange(1);
yLen = pcRange(4) - pcRange(3);
xStep = 0.3; yStep = 0.3;
nX8 = ceil(xLen / xStep / 8) * 8;
nY8 = ceil(yLen / yStep / 8) * 8;
vx = xLen / nX8;
vy = yLen / nY8;
voxelSize = [vx vy];
disp("✅ Step 5 complete: Voxel size = [" + vx + ", " + vy + "]");

%% 6) Estimate anchor boxes
disp("🔄 Step 6: Estimating anchor boxes from train data...");
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
disp("✅ Step 6 complete: Anchor boxes ready.");

%% 7) Build untrained PointPillars detector
detector = pointPillarsObjectDetector(pcRange, classNames, anchorBoxes, ...
                                     VoxelSize=voxelSize);
disp("✅ Step 7 complete: Detector initialized.");

%% 8) Set training options with dynamic checkpoint path
disp("🔄 Step 8: Setting training options with dynamic checkpoint path...");

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



% Update training options with proper validation
options = trainingOptions('adam', ...
    'Plots', "none", ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 1, ...
    'InitialLearnRate', 1e-4, ...
    'ValidationData', valPair, ...        
    'ValidationFrequency', 10, ...
    'ValidationPatience', 5, ...
    'LearnRateSchedule', "piecewise", ...
    'LearnRateDropPeriod', 20, ...
    'LearnRateDropFactor', 0.8, ...
    'CheckpointPath', checkpointDir, ...
    'CheckpointFrequency', 5, ...
    'CheckpointFrequencyUnit','epoch', ...
    'ResetInputNormalization', false, ...
    'Shuffle', "every-epoch", ...
    'BatchNormalizationStatistics','moving');

disp("✅ Step 8 complete: Training options set.");
disp("📁 Checkpoint directory: " + checkpointDir);
disp("🔍 Validation setup:");
disp("   - Training samples: " + numel(gTruthTrain) + " (with augmentation)");
disp("   - Validation samples: " + numel(gTruthVal) + " (no augmentation)");
disp("   - Max epochs: 100, early stopping after 50 epochs if no improvement");

% === Save training options to text for reference ===
optFile = fullfile(checkpointDir, 'trainingOptions.txt');
fid = fopen(optFile, 'w');
if fid > 0
    fprintf(fid, "Training options:\n\n");
    disp(options);                          % show in console
    optStr = evalc('disp(options)');        % capture as string
    fprintf(fid, "%s\n", optStr);
    
    % Also save dataset information
    fprintf(fid, "\n\nDataset Information:\n");
    fprintf(fid, "Total samples: %d\n", numel(gTruthAll));
    fprintf(fid, "Training samples: %d (%.1f%%)\n", numel(gTruthTrain), 100*numel(gTruthTrain)/numel(gTruthAll));
    fprintf(fid, "Validation samples: %d (%.1f%%)\n", numel(gTruthVal), 100*numel(gTruthVal)/numel(gTruthAll));
    fprintf(fid, "Test samples: %d (%.1f%%)\n", numel(gTruthTest), 100*numel(gTruthTest)/numel(gTruthAll));
    fprintf(fid, "\nClass names: %s\n", strjoin(classNames, ', '));
    
    fclose(fid);
    disp("📝 Saved training options to " + optFile);
else
    warning("Could not create %s", optFile);
end

%% 9) Train the model
disp("🚀 Step 9: Starting model training with augmentation...");
disp("📁 Training checkpoints will be saved to: " + checkpointDir);

[detector, info] = trainPointPillarsObjectDetector(cds, detector, options);
disp("✅ Step 9 complete: Model training finished.");

%% 10) Save model with timestamp
modelFilename = sprintf('trainedPointPillarsDetector_%s.mat', timestamp);
modelPath = fullfile(rootDir, modelFilename);

% Create training metadata
trainingMetadata = struct();
trainingMetadata.timestamp = timestamp;
trainingMetadata.checkpointDir = checkpointDir;
trainingMetadata.trainSamples = numel(gTruthTrain);
trainingMetadata.valSamples = numel(gTruthVal);
trainingMetadata.testSamples = numel(gTruthTest);
trainingMetadata.classNames = classNames;
trainingMetadata.pcRange = pcRange;
trainingMetadata.voxelSize = voxelSize;

save(modelPath, "detector", "info", "trainingMetadata");
disp("💾 Model saved to: " + modelPath);
disp("📁 Training checkpoints available at: " + checkpointDir);

% Display final training statistics if available
if ~isempty(info) && isfield(info, 'TrainingLoss')
    disp("📊 Training completed:");
    disp("   - Final training loss: " + info.TrainingLoss(end));
    if isfield(info, 'ValidationLoss') && ~isempty(info.ValidationLoss)
        disp("   - Final validation loss: " + info.ValidationLoss(end));
        disp("   - Best validation loss: " + min(info.ValidationLoss));
    end
    disp("   - Epochs completed: " + length(info.TrainingLoss));
end

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

function rowOut = sanitizeRow(rowIn)
% rowIn/rowOut = {ptCloud, boxes, labels}
pc = rowIn{1}; boxes = rowIn{2}; labels = rowIn{3};

if isempty(boxes)
    rowOut = rowIn; 
    return
end

% Remove NaN/Inf
bad = any(~isfinite(boxes), 2);

% Enforce strictly positive [L W H] if present (cols 4:6)
if size(boxes,2) >= 6
    bad = bad | boxes(:,4)<=0 | boxes(:,5)<=0 | boxes(:,6)<=0;
end

% Drop bad ones
boxes(bad,:) = [];
if ~isempty(labels)
    labels(bad) = [];
end

rowOut = {pc, boxes, labels};
end