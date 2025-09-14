% updateLabelPaths.m
% Fix paths for multiple groundTruthLidar .mat files and print updated paths.

labelDir = "/Volumes/T7_Shield/WildPose/WildPose_Proj/label/relativePath";

% List only *.mat files, then remove any that start with a dot
matFiles = dir(fullfile(labelDir, "*.mat"));
matFiles = matFiles(~startsWith({matFiles.name}, '.'));

proj = currentProject;
root = proj.RootFolder;

oldBase = "/Volumes/T7_Shield/WildPose/WildPose_Proj/";
newBase = fullfile(root);
altPaths = [oldBase, newBase];

for k = 1:numel(matFiles)
    matFile = fullfile(matFiles(k).folder, matFiles(k).name);
    S = load(matFile, 'gTruth');
    gTruth = S.gTruth;

    unresolved = changeFilePaths(gTruth, altPaths);

    % Warn if any paths couldn't be updated
    if ~isempty(unresolved)
        warning("  [!] Unresolved paths in %s:", matFiles(k).name);
        fprintf("    %s\n", unresolved);
    end

    % Print out all the paths of the updated data sources
    ds = gTruth.DataSource;
    fprintf("  Updated paths in %s:\n", matFiles(k).name);
    fprintf("    %s\n", ds.SourceName);

    % Save back updated gTruth
    save(matFile, '-struct', 'S', 'gTruth');
end

disp('Done updating all .mat files.');
