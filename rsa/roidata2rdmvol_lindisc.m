% given a set of ROIs (searchlight spheres or masks) and a pilab GLM
% instance, fit a linear discriminant for each ROI and compute some
% discriminant dissimilarity measure according to parameters in the ts
% struct.
%
% roidata2rdmvol is used when you have one response pattern per condition
% and want some pdist-based distance metric on the patterns. This function
% is used when you have raw data and a design matrix and want to take
% advantage of the better error covariance matrix estimated from the full
% model residual or to cross-validate a split-data RDM, e.g. for LDAt RDMs.
%
% named varargs:
%
% (these get passed to vol2glm_batch)
% sgolayK: degree of optional Savitsky-Golay detrend 
% sgolayF: filter size of optional Savitsky-Golay detrend
% covariatedeg: polynomial detrend degree (or 'adaptive')
% targetlabels: labels to explicitly include (default all)
% ignorelabels: labels to explicitly exclude (default none)
% glmclass: char defining CovGLM sub-class (e.g. 'CovGLM')
% glmvarargs: any additional arguments for GLM (e.g. k for RidgeGLM)
%
% (these are used in this function)
% split: indices to define GLM cvgroups (if crossvalidate) - NB NOT used in
%   vol2glm_batch
% sterrunits: scale by standard error (default false)
% crossvalidate: get split data discriminant (default false)
% 
% The gist of these options are:
% ~split && ~sterrunits = mahalanobis distance
% ~split && sterrunits = 'whitened euclidean' (?)
% split && ~sterrunits = mahalanobis classifier distance
% split && sterrunits = LDAt RDM
%
% disvol = roidata2rdmvol_lindisc(rois,designvol,epivol,[varargin])
function disvol = roidata2rdmvol_lindisc(rois,designvol,epivol,varargin)

ts = varargs2structfields(varargin,struct('sgolayK',[],'sgolayF',[],...
    'split',[],'covariatedeg',[],'targetlabels',{},'ignorelabels',{},...
    'glmclass','GLM','glmvarargs',{},'sterrunits',false,'crossvalidate',...
    false));

if ~iscell(ts.ignorelabels)
    ts.ignorelabels = {ts.ignorelabels};
end

if ~iscell(ts.split)
    % so we can easily deal to cvgroup field later
    ts.split = num2cell(ts.split);
end

% preallocate result (dissimilarity by roi)
% uh oh - it's actually not trivial knowing how many predictors we have
% now.
if isempty(ts.targetlabels)
    n = designvol.nfeatures - numel(ts.ignorelabels);
else
    % probably not terribly robust. 
    n = numel(ts.targetlabels) - numel(ts.ignorelabels);
    warning('uncharted waters - targetlabels defines dissimilarities size');
end

dissimilarities = NaN([nchoosek(n,2) rois.nsamples]);

% check for nans
nanmask = ~any(isnan(epivol.data),1);

% compute result
% to make this a parfor - need to put all the parameter setting etc outside
% the loop. Ideally we should just call a function / method on every
% iteration.
parfor n = 1:rois.nsamples
    % skip empty rois (these come out as NaN)
    validvox = full(rois.data(n,:)~=0) & nanmask;
    if ~any(validvox)
        % empty roi
        continue
    end

    % make the GLM instance for this ROI - always a cell array with one
    % entry (since we disable the split parameter under the assumption that
    % you've already used it to do a splitvol before this)
    % we skip nancell output because the data are assumed to be nan-less if
    % they made it past the validvox check
    model = vol2glm_batch(designvol,epivol(:,validvox),...
        'sgolayK',ts.sgolayK,'sgolayF',ts.sgolayF,'split',[],...
        'covariatedeg',ts.covariatedeg,'targetlabels',ts.targetlabels,...
        'ignorelabels',ts.ignorelabels,'glmclass',ts.glmclass,...
        'glmvarargs',ts.glmvarargs);
    assert(numel(model)==1,'multiple glm instances')
    % nb model can still have multiple array entries
    model = model{1};

    % fit the discriminants 
    % pairwise contrasts
    conmat = allpairwisecontrasts(model(1).npredictors);
    % probably all this if checking should be done to set up some sort of
    % constructor outside the loop. A very simple class or even a handle to
    % on of multiple sub-functions would do it.
    if ts.sterrunits
        testmeth = 'infotmap';
    else
        testmeth = 'infomahalanobis';
    end
    if ts.crossvalidate
        % split defines crossvalidation split in GLM (NB in other contexts
        % split may get passed to vol2glm_batch instead to make one GLM
        % instance per split).
        [model.cvgroup] = ts.split{:};
        cvres = cvclassificationrun(model,'discriminant',testmeth,[],conmat);
        % result - mean across splits
        res = mean(cvres,3);
    else
        % just self-fit
        w = discriminant(model,conmat);
        res = feval(testmeth,model,w,conmat);
    end
    dissimilarities(:,n) = res;
end

% convert to volume - here it is a problem that the result may have
% different nfeatures than the mask (e.g. for ROI analysis or when we do
% not run all possible searchlight spheres)
if rois.nsamples == rois.nfeatures
    % simple case - assume that samples and features are in register
    disvol = MriVolume(dissimilarities,rois,'metafeatures',struct(...
        'names',{rois.meta.samples.names}));
else
    % complicated case - need to forget the mask and write out a mask-less
    % volume. But save coordinates of ROIs to enable sanity checks later
    coords = cell(1,rois.nsamples);
    nvox = NaN([1 rois.nsamples]);
    for c = 1:rois.nsamples
        % compute centre of mass for this ROI
        coords{c} = round(mean(rois.linind2coord(rois.linind(...
            rois.data(c,:)~=0)),2));
        nvox(c) = sum(rois.data(c,:)~=0);
    end
    % make a mask-less volume 
    disvol = MriVolume(dissimilarities,[],'metafeatures',struct(...
        'names',{rois.meta.samples.names'},'centreofmass',{coords},...
        'nfeatures',nvox),'header',rois.header);
end
