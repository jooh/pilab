% batch convert designvol / epivol to a cell of GLM instances. Used
% internally in aamod_pilab_glmfit and aamod_pilab_rdms_ldt.
%
% ts fields:
%
% sgolayK: degree of optional Savitsky-Golay detrend 
% sgolayF: filter size of optional Savitsky-Golay detrend
% split: indices to define which runs go into which glmcell entry
% covariatedeg: polynomial detrend degree (or 'adaptive')
% targetlabels: labels to explicitly include (default all)
% ignorelabels: labels to explicitly exclude (default none)
% glmclass: char defining CovGLM sub-class (e.g. 'CovGLM')
% glmvarargs: any additional arguments for GLM (e.g. k for RidgeGLM)
%
% [glmcell,nancell] = vol2glm_batch(designvol,epivol,ts)
function [glmcell,nancell] = vol2glm_batch(designvol,epivol,ts)

nchunks = epivol.desc.samples.nunique.chunks;

% detrend before splitting
if ~isempty(ts.sgolayK)
    assert(isempty(ts.covariatedeg),['probably should not ' ...
        'detrend and include trend covariates at the same time']);
    fprintf('detrending (K=%d,F=%d)\n',...
        ts.sgolayK,...
        ts.sgolayF);
    epivol.sgdetrend(ts.sgolayK,...
        ts.sgolayF);
    designvol.sgdetrend(ts.sgolayK,...
        ts.sgolayF);
end

% find a reasonable polynomial degree
if strcmp(ts.covariatedeg,'adaptive')
    covdegs = NaN([1 nchunks]);
    for c = 1:nchunks
        nvol = sum(epivol.meta.samples.chunks == ...
            epivol.desc.samples.unique.chunks(c));
        % Kay's rule for deciding on n covariates based on run
        % duration.
        covdegs(c) = round(nvol * epivol.frameperiod / 60 / 2);
    end
    assert(~any(isnan(ts.covariatedeg)),...
        'failed to identify covariatedeg');
    % if we found a single deg that works for all runs, life is
    % easy
    if all(covdegs(1)==covdegs)
        ts.covariatedeg = covdegs(1);
        fprintf('adaptively selected covariate degree: %d\n',...
            ts.covariatedeg);
    else
        % otherwise we need to add new functionality to CovGLM to
        % support different covariatedeg for different runs.
        error(['Adaptive covariate deg selection failed. ' ...
            'Different covariates selected per run: ' ...
            mat2str(covdegs)]);
    end
end

% find correct labels
nlabels = length(designvol.desc.features.unique.labels);
if isempty(ts.targetlabels)
    coninds = 1:nlabels;
else
    coninds = findStrInArray(designvol.desc.features.unique.labels,...
        ts.targetlabels);
end

if ~isempty(ts.ignorelabels)
    ignoreinds = findStrInArray(designvol.desc.features.unique.labels,...
        ts.ignorelabels);
    coninds = setdiff(coninds,ignoreinds);
end

assert(~isempty(coninds),...
    'found no labels matching targetlabels/ignorelabels %s/%s',...
    ts.targetlabels,ts.ignorelabels);
ncon = length(coninds);
fprintf('storing estimates for %d conditions\n',ncon);

% empty cells don't get read properly
if isempty(ts.glmvarargs)
    ts.glmvarargs = {};
end

% split the volumes
[designcell,epicell] = splitvol(ts.split,designvol,epivol);

glmcell = cell(nsplit,1);
nancell = cell(nsplit,1);
for s = 1:nsplit
    fprintf('fitting split %d of %d...\n',s,nsplit);
    % split-specific data
    splitind = split==usplit(s);
    splitepi = epicell{s};
    splitdesign = designcell{s};

    % now mask out any NaN features
    nanmask = ~any(isnan(splitepi.data),1);
    % if you haven't NaNed out all volumes for a given feature
    % something is likely broken
    assert(all(all(isnan(splitepi.data(:,~nanmask)))),'inconsistent NaN mask detected');
    nnans = sum(~nanmask);
    if ~all(nanmask)
        fprintf(['removed %d NaN features from analysis ' ...
            '(%.2f%% of total).\n'],nnans,...
            100*(nnans/length(nanmask)));
        splitepi = splitepi(:,nanmask);
    end
    % implement some kind of glm
    glmcell{s} = vol2glm(splitdesign,splitepi,ts.glmclass,...
        ts.covariatedeg,ts.glmvarargs{:});
    nancell{s} = nanmask;
end
