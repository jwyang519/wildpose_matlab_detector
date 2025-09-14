function out_file = fix_label_typo(in_file, oldName, newName, out_file)
% Rename a label (e.g., 'Gemsbok' -> 'Springbok') in a groundTruthLidar export.
% Works whether you start from a .mat containing gTruth, or from a .m that you RUN first.
%
% Usage:
%   % If you have a .mat with a gTruth variable:
%   fix_label_typo('scene01_gTruth.mat','Gemsbok','Springbok');
%
%   % If you have a .m that constructs gTruth:
%   run('scene01_gTruth.m');                % this should create gTruth in workspace
%   save('scene01_gTruth_tmp.mat','gTruth') % save to mat
%   fix_label_typo('scene01_gTruth_tmp.mat','Gemsbok','Springbok','scene01_gTruth_fixed.mat');

    if nargin < 4
        [p,b,~] = fileparts(in_file);
        out_file = fullfile(p, [b '_fixed.mat']);
    end

    S = load(in_file);
    if ~isfield(S,'gTruth')
        error('File must contain variable "gTruth" (groundTruthLidar).');
    end
    g = S.gTruth;

    % Pull defs and data
    try
        labelDefs = g.LabelDefinitions;   % table from the object
    catch
        error('Could not read LabelDefinitions from gTruth.');
    end
    LD = g.LabelData;                     % table: one variable per label name

    % Normalize names
    oldName = string(oldName); newName = string(newName);

    haveOld = ismember(oldName, string(LD.Properties.VariableNames));
    haveNew = ismember(newName, string(LD.Properties.VariableNames));

    % --- 1) Update label definitions table
    defNames = string(labelDefs.Name);
    if any(defNames == oldName)
        if any(defNames == newName)
            % If both exist, drop the old row
            labelDefs(defNames == oldName, :) = [];
        else
            % Just rename the row
            labelDefs.Name(defNames == oldName) = cellstr(newName);
        end
    end

    % --- 2) Update LabelData columns
    if haveOld && haveNew
        % Merge per-frame tables: append old rows into new, keep variable union safely
        LD.(char(newName)) = cellfun(@(A,B) cat_tables(B, A), LD.(char(newName)), LD.(char(oldName)), 'UniformOutput', false);
        LD.(char(oldName)) = []; % drop the mistaken column
    elseif haveOld && ~haveNew
        % Simple rename of the column
        LD.Properties.VariableNames{ strcmp(LD.Properties.VariableNames, char(oldName)) } = char(newName);
    else
        warning('No column named "%s" found in LabelData. Nothing to rename there.', oldName);
    end

    % --- 3) Rebuild a clean groundTruthLidar and save
    g2 = groundTruthLidar(g.DataSource, labelDefs, LD);
    gTruth = g2; %#ok<NASGU>
    save(out_file,'gTruth');
    fprintf('✔ Saved fixed file: %s\n', out_file);
end

function T = cat_tables(Tkeep, Tadd)
% Append Tadd rows into Tkeep, aligning variables safely (handles empty).
    if isempty(Tadd)
        T = Tkeep; return;
    end
    if isempty(Tkeep)
        T = Tadd; return;
    end
    % Align variable sets (outer union)
    allVars = unique([Tkeep.Properties.VariableNames, Tadd.Properties.VariableNames]);
    Tkeep = ensure_vars(Tkeep, allVars);
    Tadd  = ensure_vars(Tadd , allVars);
    T = [Tkeep; Tadd];
end

function T = ensure_vars(T, allVars)
% Add any missing variables as empty columns with appropriate types.
    cur = T.Properties.VariableNames;
    miss = setdiff(allVars, cur);
    for k = 1:numel(miss)
        % Best-effort: add missing column as empty column of type double; tables will accept it.
        T.(miss{k}) = repmat(missingLike(T, miss{k}), height(T), 1);
    end
end

function v = missingLike(~, ~)
% Return a scalar missing value; using double NaN is acceptable for placeholder.
    v = NaN; % Lidar Labeler cuboid tables ignore extra NaN cols; they won't be saved back anyway.
end