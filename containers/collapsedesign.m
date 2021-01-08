% collapse regressors in designvol according to cell array levels. For each
% element in label we use strfindcell to find matching
% designvol.meta.features.labels, which are summed into a new design
% matrix. Regressors that do not match any of the search patterns are
% appended at the end of the new design matrix.
%
% newdesignvol = collapsedesign(designvol, levels)
function newdesignvol = collapsedesign(designvol, levels)

labels = designvol.meta.features.labels(:);
targetind = [];
newlabels = {};
for n = 1:numel(levels)
    targetind = [targetind strfindcell(labels,levels{n})];
    newlabels = [newlabels strrep(labels(targetind(:,end)),levels{n},'collapsed')];
end
allgood = arrayfun(@(x)isequal(newlabels{x,:}),1:size(newlabels,1),'uniformoutput',1);
assert(all(allgood), 'regressors are not sorted properly, try again');

designmatrix = designvol.data;
newdesignmatrix = zeros(size(designmatrix,1), size(targetind,1));
for regind = 1:size(newdesignmatrix,2)
    newdesignmatrix(:,regind) = sum(designmatrix(:,targetind(regind,:)),2);
end

% add on the stuff that's unaffected by collapsing
[~,newind] = setdiff(1:size(designmatrix,2),targetind);
newdesignmatrix = [newdesignmatrix, designmatrix(:,newind)];
newlabels = [newlabels(:,1); labels(newind')];

% build output volume
newdesignvol = feval(class(designvol),newdesignmatrix,'frameperiod',...
    designvol.frameperiod, 'metafeatures', struct('labels', {newlabels'}), ...
    'metasamples', designvol.meta.samples);