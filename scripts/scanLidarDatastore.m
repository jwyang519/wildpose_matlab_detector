function report = scanLidarDatastore(ds, name, classNames)
% Scan a (possibly transformed/combined) lidar datastore.
% Flags only problematic rows; prints a short summary and returns a table.

    if nargin < 3, classNames = []; end
    reset(ds);

    nRows = 0; bad = 0; rep = struct([]);   % accumulate issues

    while hasdata(ds)
        nRows = nRows + 1;
        row = read(ds);

        probs = strings(0,1);

        % Expect row = {pointCloud, boxes(Mx9), labels(Mx1 categorical)}
        okForm = iscell(row) && numel(row) == 3;
        if ~okForm
            probs(end+1) = "row not 1x3 cell";
            bx = []; lb = [];
        else
            bx = row{2};
            lb = row{3};
        end

        % --- Boxes checks ---
        if ~isempty(bx)
            if ~(ismatrix(bx) && size(bx,2) == 9)
                probs(end+1) = sprintf("boxes not Mx9 (got %s)", mat2str(size(bx)));
            else
                if any(~isfinite(bx(:)))
                    probs(end+1) = "boxes contain NaN/Inf";
                end
                if any(bx(:,4:6) <= 0, 'all')
                    probs(end+1) = "nonpositive xlen/ylen/zlen";
                end
            end
        end

        % --- Labels checks ---
        if ~isempty(bx)  % when there are boxes, labels count must match M
            if numel(lb) ~= size(bx,1)
                probs(end+1) = sprintf("labels length (%d) != num boxes (%d)", numel(lb), size(bx,1));
            end
        end
        if ~isempty(lb)
            if ~iscategorical(lb)
                probs(end+1) = "labels not categorical";
            elseif ~isempty(classNames)
                % ensure membership in your declared class set
                ok = ismember(lb, categorical(string(classNames)));
                if any(~ok)
                    badCats = unique(string(lb(~ok)));
                    probs(end+1) = "label(s) outside classNames: " + strjoin(badCats, ", ");
                end
            end
        end

        % --- record any problems ---
        if ~isempty(probs)
            bad = bad + 1;
            rep(end+1).rowIndex = nRows; %#ok<AGROW>
            rep(end).boxSize    = size(bx);
            rep(end).nLabels    = numel(lb);
            rep(end).problems   = strjoin(probs, " | ");
        end
    end

    fprintf('[%s] scanned %d rows. bad=%d (%.1f%%)\n', name, nRows, bad, 100*bad/max(1,nRows));

    if ~isempty(rep)
        report = struct2table(rep);
        disp(report(1:min(20,height(report)), :));   % show first few bad rows
    else
        report = table();
        disp("No issues found ✅");
    end
end