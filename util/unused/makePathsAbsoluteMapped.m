function makePathsAbsoluteMapped(labelDir, rootPairs, pattern, overwrite)
% Map stored DataSource paths to absolute, existing paths and save.
% labelDir  : folder with .mat (groundTruthLidar)
% rootPairs : N×2 string/char array [storedRoot, actualRootOnThisOS]
% pattern   : optional, default "*.mat"
% overwrite : optional, default false (writes *_abs.mat)

if nargin < 3 || isempty(pattern),  pattern = "*.mat"; end
if nargin < 4, overwrite = false; end

labelDir  = string(labelDir);
rootPairs = string(rootPairs);

L = dir(fullfile(labelDir, pattern));
mask = ~startsWith({L.name},'.') & ~startsWith({L.name},'._') & ~strcmp({L.name},'.DS_Store');
L = L(mask);

fprintf("Scanning %d MAT files under %s\n", numel(L), labelDir);

for k = 1:numel(L)
    inFile = fullfile(L(k).folder, L(k).name);
    S = load(inFile);  % preserve original var names

    [g,gvar] = findGTL(S);
    if isempty(g)
        fprintf("[SKIP] %s (no groundTruthLidar)\n", L(k).name);
        continue;
    end

    srcs = string(g.DataSource.SourceName);
    if isempty(srcs)
        fprintf("[SKIP] %s (empty SourceName)\n", L(k).name);
        continue;
    end

    % Build N×2 CELL of char: {p_current, p_new_abs}
    altPairs = cell(0,2);
    changed  = 0;

    for r = 1:numel(srcs)
        pcur = srcs(r);                         % exact stored value
        pcurNorm = normpath_str(pcur);

        newAbs = "";
        for p = 1:size(rootPairs,1)
            fromRoot = normpath_str(rootPairs(p,1));
            toRoot   = normpath_str(rootPairs(p,2));
            if startsWith(pcurNorm, fromRoot)
                newAbs = replace(pcurNorm, fromRoot, toRoot);
                break;
            end
        end

        % If no mapping matched but the stored path exists as-is, keep it.
        if newAbs == "" && isfolder(char(pcur))
            newAbs = pcurNorm;
        end

        if newAbs ~= "" && isfolder(char(newAbs))
            altPairs(end+1,1) = {char(pcur)};      %#ok<AGROW>
            altPairs(end  ,2) = {char(newAbs)};
            changed = changed + 1;
        else
            fprintf("  [WARN] %s: cannot resolve to existing folder:\n         %s\n", ...
                L(k).name, pcur);
            if newAbs ~= ""
                fprintf("         (mapped to) %s  [NOT FOUND]\n", newAbs);
            end
        end
    end

    if isempty(altPairs)
        fprintf("[SKIP] %s (no resolvable mappings)\n", L(k).name);
        continue;
    end

    unresolved = changeFilePaths(g, altPairs);   % updates in-place
    fprintf("[OK] %s: requested %d changes, unresolved %d\n", ...
        L(k).name, changed, numel(unresolved));
    if ~isempty(unresolved)
        disp("  e.g. unresolved:"); disp(unresolved(1:min(3,end)));
    end

    % Save back
    S.(gvar) = g;
    if overwrite
        save(inFile, "-struct", "S", "-v7.3");
        fprintf("      Saved IN-PLACE.\n");
    else
        [p,b,e] = fileparts(inFile);
        outFile = fullfile(p, b + "_abs" + e);
        save(outFile, "-struct", "S", "-v7.3");
        fprintf("      Wrote %s\n", outFile);
    end
end
end

% --- helpers ---
function [g,name] = findGTL(S)
g = []; name = '';
for f = fieldnames(S).'
    if isa(S.(f{1}), 'groundTruthLidar'), g=S.(f{1}); name=f{1}; return; end
end
end

function s = normpath_str(s)
s = string(s);
s = replace(s, '\', '/');  % normalize for matching
s = regexprep(s, '/+$','');
end