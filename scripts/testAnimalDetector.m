%% 1) Detect on test set – keep one row per frame
numF          = height(testData);

detBoxesCell  = cell(numF,1);
detScoresCell = cell(numF,1);
detLabelsCell = cell(numF,1);

for k = 1:numF
    ptCloud = pcread(testData.ptCloudFilename{k});

    % Run detector (full 3-D cuboids -> pick BEV 5-D rotated rects)
    [bboxes,scores,labels] = detect(detector,ptCloud,Threshold=0.01);

    detBoxesCell{k}  = bboxes(:,[1 2 4 5 9]);   % [xc yc w h yaw]
    detScoresCell{k} = scores;
    detLabelsCell{k} = labels;                  % already categorical
end

detectionResults = table(detBoxesCell,detScoresCell,detLabelsCell, ...
    'VariableNames',{'Boxes','Scores','Labels'});

%% 2) Build matching ground-truth table (same row order)
gtBoxesCell   = cell(numF,1);
gtLabelsCell  = cell(numF,1);

for k = 1:numF
    cub   = testData{k,2}{1};          % original 9-D boxes
    gtBoxesCell{k}  = cub(:,[1 2 4 5 9]);       % BEV projection
    nB = size(cub,1);
    gtLabelsCell{k} = repmat(categorical("Giraffe"), nB, 1);
end

groundTruthData = table(gtBoxesCell,gtLabelsCell, ...
    'VariableNames',{'Boxes','Labels'});

%% 3) Evaluate (empty rows are OK)
metrics = evaluateObjectDetection(detectionResults,groundTruthData, ...
    AdditionalMetrics="AOS");

[dsSummary, classSummary] = summarize(metrics)