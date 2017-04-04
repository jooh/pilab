% convert between two internal formats used by pilab roidata analysis functions:
% subres is a struct array with one field per subject (e.g., output from
% roidata_rsa), groupres is a scalar struct with subject data stacked in the
% third dimension. For usage, see e.g.  roidata_rfx.
%
% subres: struct array input with results from e.g. roidata_rsa, roidata_glm. If
%   you want to code your own, see these functions for conventions
% targetfield: ('r') field with relevant data in subres
% assumeregister: (false) assume that all rois are present in the same order in
%   all subres entries. Speeds things up, and preserves ordering of ROIs.
%
% groupres = subres2groupres(subres,targetfield='r',assumeregister=false)
function groupres = subres2groupres(subres,targetfield,assumeregister)

% input check
if ~exist('targetfield','var') || isempty(targetfield)
    targetfield = 'r';
end
assert(isfield(subres,targetfield),'targetfield %s is missing in subres',...
    targetfield);
if ~exist('assumeregister','var') || isempty(assumeregister)
    assumeregister = false;
end

if numel(subres)==1 && isfield(subres,'z_subject')
    % already groupres
    groupres = subres;
    return
end

targets = intersect(fieldnames(subres),...
    {targetfield,'pperm','ppara','medianboot','medianste',...
    'nfeatures','tail'});

% we assume that contrasts are identical.  we could support
% this case in theory but let's not for now since it likely
% indicates a user error.
concell = {subres.rows_contrast};
assert(isequal(concell{1},concell{:}),...
    'different predictors in different subjects');
ncon = length(concell{1});

% just take rois from first subject
urois = subres(1).cols_roi;
if ~assumeregister
    % find all possible ROIs (not all subjects will have all ROIs
    % necessarily)
    allrois = {subres.cols_roi};
    allrois = horzcat(allrois{:});
    allurois = unique(allrois);
    % sort in same order as subject 1's ROIs as far as possible.
    urois = [urois setdiff(allurois,urois)];
end
nroi = length(urois);
nsub = length(subres);
dat = NaN([ncon nroi nsub]);
groupres = struct('rows_contrast',{subres(1).rows_contrast},...
    'cols_roi',{urois},targetfield,dat,'nfeatures',NaN([1 nroi nsub]),...
    'tail',{subres(1).tail},'z_subject',{{subres.name}});

% populate the groupres struct
if assumeregister
    % hey, this is easy
    groupres = collapsestruct(subres,@zcat);
    assert(~any(isnan(groupres.(targetfield)(:))),...
        'nans in targetfield not supported in assumeregister mode');
else
    % have to be careful not to average over different rois
    for s = 1:nsub
        for r = 1:length(subres(s).cols_roi)
            thisroi = subres(s).cols_roi{r};
            indroi = strcmp(groupres.cols_roi,thisroi);
            for t = targets(:)'
                tstr = t{1};
                if strcmp(tstr,'tail')
                    assert(isequal(groupres.tail,subres.tail),...
                        'mismatched tails across subjects');
                else
                    groupres.(tstr)(:,indroi,s) = ...
                        subres(s).(tstr)(:,r);
                end
            end
        end
    end
end

if isfield(subres,'custom')
    ncustom = numel(subres(1).custom);
    groupres.custom = cell(size(subres(1).custom));
    for s = 1:nsub
        assert(isequal(ncustom,numel(subres(s).custom)),...
            'mismatched numel(custom) across subjects');
        % work out the ROI indices for this subject
        [~,ginds,sinds] = intersect(groupres.cols_roi,subres(s).cols_roi,...
            'stable');
        for c = 1:ncustom
            assert(ismatrix(subres(s).custom{c}),['>2d custom fields '...
                'are not supported']);
            groupres.custom{c} = zcat(groupres.custom{c},...
                NaN([size(subres(s).custom{c},1),nroi]));
            % so we make no assumptions about the shape of custom here except
            % that ROIs are in second dim
            groupres.custom{c}(:,ginds,end) = subres(s).custom{c}(:,sinds);
        end
    end
end
