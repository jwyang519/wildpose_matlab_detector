function S = summarize_gtruth(matFile)
% Quick summary of labels in a groundTruthLidar MAT file.
% Usage: summarize_gtruth('scene01_gTruth.mat')

G = load(matFile,'gTruth'); g = G.gTruth;
names = string(g.LabelDefinitions.Name);
types = string(g.LabelDefinitions.Type);

% Count instances per label (sum of rows in each per-frame table)
counts = zeros(numel(names),1);
for i = 1:numel(names)
    col = char(names(i));
    if ismember(col, g.LabelData.Properties.VariableNames)
        counts(i) = sum(cellfun(@(t) ~isempty(t) * height(t), g.LabelData.(col)));
    end
end

% Frames with any label
nFrames = height(g.LabelData);
hasAny = false(nFrames,1);
for i = 1:numel(names)
    col = char(names(i));
    if ismember(col, g.LabelData.Properties.VariableNames)
        cells = g.LabelData.(col);
        hasAny = hasAny | cellfun(@(t) ~isempty(t) && height(t) > 0, cells);
    end
end

% (Optional) list attributes defined per label (if present)
attrStr = strings(numel(names),1);
if ismember("Attributes", g.LabelDefinitions.Properties.VariableNames)
    for i = 1:numel(names)
        A = g.LabelDefinitions.Attributes{i};
        if istable(A) && ~isempty(A) && ismember("Name", A.Properties.VariableNames)
            attrStr(i) = strjoin(string(A.Name), ',');
        else
            attrStr(i) = "";
        end
    end
end

T = table(names, types, counts, attrStr, ...
    'VariableNames', {'Label','Type','InstanceCount','Attributes'});

fprintf('File: %s\n', matFile);
fprintf('Frames: %d | Frames with any label: %d\n', nFrames, nnz(hasAny));
fprintf('Total instances: %d\n\n', sum(counts));
disp(T);

% Return struct if you want to reuse programmatically
S = struct('Frames', nFrames, 'FramesWithLabels', nnz(hasAny), ...
           'TotalInstances', sum(counts), 'PerLabel', T);
end