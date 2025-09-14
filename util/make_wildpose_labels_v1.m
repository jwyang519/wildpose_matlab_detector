% make_wildpose_labels_master.m
function out = make_wildpose_labels_master()
    species = {'Cheetah','Giraffe','Ostrich','Gemsbok','Springbok'}; % start set

    ldc = labelDefinitionCreatorLidar;         % <— Lidar Labeler creator (no SignalType column)
    for k = 1:numel(species)
        addLabel(ldc, species{k}, labelType.Cuboid, 'Group','Animal');
    end
    % Example ROI attribute (optional)
    % addAttribute(ldc,'Cheetah','Occluded',attributeType.Logical,false);

    labelDefs = create(ldc);
    save('wildpose_labels_master.mat','labelDefs');   % <- import this in the app
    out = 'wildpose_labels_master.mat';
    fprintf('Wrote %s\n', out);
end