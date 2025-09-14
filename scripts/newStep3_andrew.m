%% 3) Stratified split into train/val/test from fullData table
% disp("🔄 Step 3: Spliting data 70/30...");
% rng(1);
% shuffledIndices = randperm(size(processedLabels,1));
% idx = floor(0.7 * length(shuffledIndices));
% 
% trainData = croppedPointCloudObj(shuffledIndices(1:idx),:);
% testData = croppedPointCloudObj(shuffledIndices(idx+1:end),:);
% 
% trainLabels = processedLabels(shuffledIndices(1:idx),:);
% testLabels = processedLabels(shuffledIndices(idx+1:end),:);
% 
% writeFiles = true;
dataLocation = fullfile(outputFolder,'InputData');
% [trainData,trainLabels] = saveptCldToPCD(trainData,trainLabels,...
%     dataLocation,writeFiles);
load("trainlabels.mat", "trainLabels");
load("classNames.mat", "classNames");
% lds = fileDatastore(dataLocation,'ReadFcn',@(x) pcread(x));
% bds = boxLabelDatastore(trainLabels);
% cds = combine(lds,bds);
% disp("✅ Step 3: Spliting done, file datastore saved.");

Idx = [];
for i = 1:size(trainLabels,1)
    if (isempty(trainLabels.Gemsbok{i}) && isempty(trainLabels.Giraffe{i}) && isempty(trainLabels.Lion{i}) && isempty(trainLabels.Wildbeest{i}))
        Idx = [Idx;i];
    end
end
lds = fileDatastore(dataLocation,'ReadFcn',@(x) pcread(x));
fNames = lds.Files;
fNames(Idx) = [];
lds = fileDatastore(fNames,'ReadFcn',@(x) pcread(x));
trainLabels(Idx,:) = [];
bds = boxLabelDatastore(trainLabels);

cds = combine(lds,bds);
disp("✅ Step 3: Spliting done, file datastore saved.");