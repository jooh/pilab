% batch processing of session-specific RDM modelling. A wrapper around
% roidata2rdmvol for handling the case where you want to compute RDMs
% separately for different splits of the data.
%
% [disvol,splitdisvolcell] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)
function [disvol,splitdisvolcell] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)

getArgs(varargin,{'split',[],'glmvarargs',{},'cvsplit',[],...
    'glmclass',[],'sterrunits',[],'crossvalidate',[],...
    'minvoxeln',1,'batchsize',5000});

% split the data into appropriately pre-processed cell arrays
% (now skip preprocessing - we assume this has already happened)
[designcell,epicell] = splitvol(split,designvol,epivol);

% track NaN features across splits may appear in different runs if
% nan masking - nans should only really happen here if you enter
% empty ROIs
nsplit = length(designcell);
splitdiscvolcell = cell(nsplit,1);

if ~iscell(glmvarargs)
    if isempty(glmvarargs)
        glmvarargs = {};
    else
        glmvarargs = {glmvarargs};
    end
end

if ischar(cvsplit)
    cvsplit = eval(cvsplit);
end

% track NaN features - may appear in different runs if nan masking
nanmask = false([nsplit rois.nsamples]);
sumdata = [];

% run the beast
for sp = 1:nsplit
    fprintf('running rois for split %d of %d...\n',sp,nsplit);
    % cart off to new function
    splitdisvolcell{sp} = roidata2rdmvol_lindisc(rois,...
        designcell{sp},epicell{sp},...
        'split',cvsplit,...
        'glmclass',glmclass,...
        'glmvarargs',glmvarargs,'sterrunits',sterrunits,...
        'crossvalidate',crossvalidate,'minvoxeln',...
        minvoxeln,'batchsize',batchsize);
    if isempty(sumdata)
        sumdata = splitdisvolcell{sp}.data;
    else
        sumdata = sumdata + splitdisvolcell{sp}.data;
    end
    nanmask(sp,:) = any(isnan(splitdisvolcell{sp}.data),1);
end % sp 1:nsplit
% remove any nan features from all sessdisvols
anynan = any(nanmask,1);
if any(anynan)
    nnans = sum(anynan);
    fprintf(['removed %d NaN ROIs from analysis ' ...
        '(%.2f%% of total).\n'],nnans,...
        100*(nnans/length(anynan)));
end
splitdisvolcell = cellfun(@(dv)dv(:,~anynan),splitdisvolcell,...
    'uniformoutput',false);
% and from sums
sumdata(:,anynan) = [];

% extract meta features for mean rdm vol (needs to be after main
% loop to avoid nan ROIs) and write out
mfeatures = splitdisvolcell{1}.meta.features;
for sp = 1:nsplit
    if isfield(splitdisvolcell{sp}.meta.features,'nfeatures')
        fn = sprintf('nfeatures_split%02d',sp);
        mfeatures.(fn) = ...
            splitdisvolcell{sp}.meta.features.nfeatures;
    end
    if isfield(splitdisvolcell{sp}.meta.features,'centreofmass')
        fn = sprintf('centreofmass_split%02d',sp);
        mfeatures.(fn) = ...
            splitdisvolcell{sp}.meta.features.centreofmass;
    end
end


% make average RDM across sessions 
if isa(splitdisvolcell{1},'MriVolume')
    disvol = MriVolume(sumdata/nsplit,splitdisvolcell{1},...
        'metafeatures',mfeatures);
else
    disvol = BaseVolume(sumdata/nsplit,...
        'metafeatures',mfeatures);
end
