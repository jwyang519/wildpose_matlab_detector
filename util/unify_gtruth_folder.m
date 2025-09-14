function unify_gtruth_folder(rootDir, masterMat, varargin)
% Unify all groundTruthLidar files under rootDir to a single master label dictionary.
% See usage examples at bottom.

% ---------- parse options ----------
p = inputParser;
p.addRequired('rootDir', @(s)ischar(s)||isstring(s));
p.addRequired('masterMat', @(s)ischar(s)||isstring(s));
p.addParameter('RenameMap', containers.Map('KeyType','char','ValueType','char'));
p.addParameter('DropExtraLabels', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('StrictAttributes', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('Overwrite', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('Suffix', '_uniform', @(s)ischar(s)||isstring(s));
p.addParameter('DryRun', false, @(x)islogical(x)&&isscalar(x));
p.parse(rootDir, masterMat, varargin{:});
opt = p.Results;
rootDir   = char(rootDir);
masterMat = char(masterMat);
suffix    = char(opt.Suffix);

% ---------- load master labelDefs ----------
M = load(masterMat);
labelDefs = [];
for k = fieldnames(M).'
    T = M.(k{1});
    if istable(T) && ismember('Name', T.Properties.VariableNames)
        labelDefs = T; break;
    end
end
if isempty(labelDefs)
    error('Could not find a labelDefs table in %s', masterMat);
end
masterNames = string(labelDefs.Name);

% helper: attributes spec for a label
    function A = masterAttrFor(lbl)
        A = [];
        if ismember('Attributes', labelDefs.Properties.VariableNames)
            r = find(masterNames == string(lbl), 1);
            if ~isempty(r)
                A = labelDefs.Attributes{r};
                if ~isempty(A) && ~istable(A), A = []; end
            end
        end
    end

% ---------- enumerate MAT files, ignore hidden/AppleDouble ----------
files = dir(fullfile(rootDir, '**', '*.mat'));
files = files(~[files.isdir]);
names = {files.name};
files = files(~startsWith(names, '.'));   % skips ._foo.mat, etc.

fprintf('Scanning %d MAT files under %s\n', numel(files), rootDir);

hasRenameVars = exist('renamevars','file') == 2;
nChanged = 0; nSkipped = 0;

for i = 1:numel(files)
    inFile = fullfile(files(i).folder, files(i).name);
    try
        S = load(inFile);
    catch ME
        warning('Skipping unreadable MAT: %s (%s)', inFile, ME.message);
        nSkipped = nSkipped + 1; continue;
    end

    % find a groundTruthLidar object
    g = [];
    for k = fieldnames(S).'
        if isa(S.(k{1}),'groundTruthLidar'), g = S.(k{1}); break; end
    end
    if isempty(g)
        if isfield(S,'gTruth') && isa(S.gTruth,'groundTruthLidar')
            g = S.gTruth;
        else
            nSkipped = nSkipped + 1; continue;
        end
    end

    % Always work on a COPY of the label table
    LD = g.LabelData;
    haveNames = string(LD.Properties.VariableNames);
    origHave  = haveNames;

    % --- step 1: rename to canonical names via map (if provided) ---
    if ~isempty(opt.RenameMap)
        keys = string(opt.RenameMap.keys);
        for kk = 1:numel(keys)
            old = keys(kk);
            new = string(opt.RenameMap(char(old)));
            if any(haveNames == old) && old ~= new
                if hasRenameVars
                    LD = renamevars(LD, char(old), char(new));
                else
                    % fallback pre-R2020a
                    LD.(char(new)) = LD.(char(old));
                    LD.(char(old)) = [];
                end
            end
        end
        haveNames = string(LD.Properties.VariableNames);
    end

    % --- step 2: add missing labels (empty cells) ---
    missing = setdiff(masterNames, haveNames);
    for nm = missing(:).'
        LD.(char(nm)) = cell(height(LD),1);
    end

    % --- step 3: drop extras not in master (optional) ---
    extras = setdiff(haveNames, masterNames);
    if ~isempty(extras) && opt.DropExtraLabels
        LD(:, cellstr(extras)) = [];
    end

    % --- step 4: backfill/prune attributes per label column ---
    haveNames2 = string(LD.Properties.VariableNames);
    for nm = haveNames2(:).'
        wantA = masterAttrFor(nm);
        col = LD.(char(nm));
        if ~iscell(col), continue; end
        for r = 1:numel(col)
            T = col{r};
            if isempty(T) || ~istable(T), continue; end

            % add any missing attributes with defaults
            if ~isempty(wantA)
                for ak = 1:height(wantA)
                    aName = string(wantA.Name(ak));
                    if ~ismember(aName, string(T.Properties.VariableNames))
                        def = [];
                        if ismember('DefaultValue', wantA.Properties.VariableNames)
                            def = wantA.DefaultValue{ak};
                        end
                        % fill with the right shape/type
                        if islogical(def)
                            T.(char(aName)) = repmat(def, height(T), 1);
                        elseif isstring(def) || ischar(def)
                            T.(char(aName)) = repmat(string(def), height(T), 1);
                        else
                            T.(char(aName)) = repmat({def}, height(T), 1);
                        end
                    end
                end
            end

            % remove attributes not in master (if strict)
            if opt.StrictAttributes && ~isempty(wantA)
                keep = ["Position"; string(wantA.Name)];
                T = T(:, intersect(keep, string(T.Properties.VariableNames)));
            end
            col{r} = T;
        end
        LD.(char(nm)) = col;
    end

    % --- step 5: rebuild a NEW groundTruthLidar with unified defs + LD ---
    g2 = groundTruthLidar(g.DataSource, labelDefs, LD);

    % --- step 6: save ---
    if opt.DryRun
        fprintf('[DRY] %s: %d→%d labels (missing:%d, extras:%d)\n', ...
            files(i).name, numel(origHave), width(g2.LabelData), numel(missing), numel(extras));
    else
        if opt.Overwrite
            backup = [inFile '.bak'];
            if ~exist(backup,'file'), copyfile(inFile, backup); end
            gTruth = g2; %#ok<NASGU>
            save(inFile, 'gTruth');
            fprintf('[OK]   %s (overwritten). backup: %s\n', files(i).name, backup);
        else
            [p,b,~] = fileparts(inFile);
            outFile = fullfile(p, [b suffix '.mat']);
            gTruth = g2; %#ok<NASGU>
            save(outFile, 'gTruth');
            fprintf('[OK]   %s -> %s\n', files(i).name, [b suffix '.mat']);
        end
    end

    nChanged = nChanged + 1;
end

fprintf('Done. Changed: %d, skipped (no/invalid gTruth): %d\n', nChanged, nSkipped);
end