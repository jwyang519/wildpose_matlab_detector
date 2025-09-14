function makePathsRelative_smart(labelDir, projectRoot, rootPairs, pattern, overwrite)
% Convert groundTruthLidar DataSource paths to RELATIVE, robust across Win/Mac roots.
% labelDir    : folder with .mat files
% projectRoot : absolute project root (the folder you'll cd into when using the data)
% rootPairs   : N×2 string array mapping stored roots -> THIS-OS absolute roots
% pattern     : (optional) default "*.mat"
% overwrite   : (optional) default false (writes *_rel.mat alongside)

if nargin < 4 || isempty(pattern),  pattern = "*.mat"; end
if nargin < 5, overwrite = false; end

labelDir    = string(labelDir);
projectRoot = normstr(projectRoot);
rootPairs   = string(rootPairs);

L = dir(fullfile(labelDir, pattern));
mask = ~startsWith({L.name},'.') & ~startsWith({L.name},'._') & ~strcmp({L.name},'.DS_Store');
L = L(mask);

fprintf("Scanning %d MAT files under %s\n", numel(L), labelDir);

for k = 1:numel(L)
    inFile = fullfile(L(k).folder, L(k).name);
    S = load(inFile);       % keep original var names

    [g,gvar] = findGTL(S);
    if isempty(g)
        fprintf("[SKIP] %s (no groundTruthLidar)\n", L(k).name); 
        continue;
    end

    src = string(g.DataSource.SourceName);
    if isempty(src)
        fprintf("[SKIP] %s (empty SourceName)\n", L(k).name);
        continue;
    end

    % Build an N×2 string array for changeFilePaths
    ALTS = strings(0,2);
    for r = 1:numel(src)
        sExact = src(r);                      % EXACT stored value
        sNorm  = normstr(sExact);

        % Case 1: already under projectRoot -> relative directly
        if startsWith(sNorm, projectRoot)
            rel = relpath_local(sNorm, projectRoot);
            ALTS(end+1, :) = [sExact, rel]; %#ok<AGROW>
            continue;
        end

        % Case 2: translate via rootPairs, then make relative
        mapped = false;
        for p = 1:size(rootPairs,1)
            fromP = normstr(rootPairs(p,1));
            toP   = normstr(rootPairs(p,2));  % absolute on THIS OS
            if startsWith(sNorm, fromP)
                absHere = replace(sNorm, fromP, toP);
                rel     = relpath_local(absHere, projectRoot);
                ALTS(end+1, :) = [sExact, rel]; %#ok<AGROW>
                mapped = true;
                break;
            end
        end

        if ~mapped
            fprintf("  [WARN] %s: no root mapping for\n         %s\n", L(k).name, sExact);
        end
    end

    if isempty(ALTS)
        fprintf("[SKIP] %s (no mappable paths)\n", L(k).name);
        continue;
    end

    unresolved = changeFilePaths(g, ALTS);
    fprintf("[OK] %s: requested %d changes, unresolved %d\n", L(k).name, size(ALTS,1), numel(unresolved));
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
        outFile = fullfile(p, b + "_rel" + e);
        save(outFile, "-struct", "S", "-v7.3");
        fprintf("      Wrote %s\n", outFile);
    end
end
end

% --- helpers ---
function s = normstr(s)
s = replace(string(s), '\','/');      % normalize separators to '/'
s = regexprep(s, '/+$', '');          % strip trailing slashes
end

function [g,name] = findGTL(S)
g = []; name = '';
for f = fieldnames(S).'
    if isa(S.(f{1}), 'groundTruthLidar'), g=S.(f{1}); name=f{1}; return; end
end
end

function rel = relpath_local(targetAbs, baseAbs)
% Compute target path RELATIVE to base path. Both inputs are absolute strings.
t = normstr(targetAbs);
b = normstr(baseAbs);
tParts = split(t, '/');  bParts = split(b, '/');

% Find common prefix length
n = min(numel(tParts), numel(bParts));
i = 1;
while i <= n && tParts(i) == bParts(i)
    i = i + 1;
end
commonLen = i - 1;

% How many ".." to go up from base
upCount = numel(bParts) - commonLen;
if upCount > 0
    up = strjoin(repmat("..",1,upCount), "/");
else
    up = "";
end

% Down segments from target
down = strjoin(tParts(commonLen+1:end), "/");

if up == ""
    rel = down;
elseif down == ""
    rel = up;
else
    rel = up + "/" + down;
end
if rel == ""
    rel = ".";   % same folder
end
rel = string(rel);
end