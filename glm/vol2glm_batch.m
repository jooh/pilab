% batch convert designvol / epivol to a cell of GLM instances. Used
% internally in aamod_pilab_glmfit and aamod_pilab_rdms_ldt.
%
% named varargs:
%
% sgolayK: degree of optional Savitsky-Golay detrend 
% sgolayF: filter size of optional Savitsky-Golay detrend
% split: indices to define which runs go into which glmcell entry
% covariatedeg: polynomial detrend degree (or 'adaptive')
% targetlabels: labels to explicitly include (default all)
% ignorelabels: labels to explicitly exclude (default none)
% glmclass: char defining CovGLM sub-class (e.g. 'CovGLM')
% glmvarargs: any additional arguments for GLM (e.g. k for RidgeGLM)
% usegpu: initialise data as gpuArray (default false)
%
% [glmcell,nancell,predictornames] = vol2glm_batch(designvol,epivol,varargin)
function [glmcell,nancell,predictornames] = vol2glm_batch(designvol,epivol,varargin)

getArgs(varargin,{'sgolayK',[],'sgolayF',[],'split',[],...
    'covariatedeg',[],'targetlabels',{},'ignorelabels',{},...
    'glmclass','CovGLM','glmvarargs',{},'usegpu',false});

nchunks = epivol.desc.samples.nunique.chunks;

% detrend before splitting
if ~isempty(sgolayK)
    assert(isempty(covariatedeg),['probably should not ' ...
        'detrend and include trend covariates at the same time']);
    fprintf('detrending (K=%d,F=%d)\n',sgolayK,sgolayF);
    epivol.sgdetrend(sgolayK,sgolayF);
    designvol.sgdetrend(sgolayK,sgolayF);
end

% find a reasonable polynomial degree
if strcmp(covariatedeg,'adaptive')
    covdegs = NaN([1 nchunks]);
    for c = 1:nchunks
        nvol = sum(epivol.meta.samples.chunks == ...
            epivol.desc.samples.unique.chunks(c));
        % Kay's rule for deciding on n covariates based on run
        % duration.
        covdegs(c) = round(nvol * epivol.frameperiod / 60 / 2);
    end
    assert(~any(isnan(covariatedeg)),...
        'failed to identify covariatedeg');
    % if we found a single deg that works for all runs, life is
    % easy
    if all(covdegs(1)==covdegs)
        covariatedeg = covdegs(1);
        %fprintf('adaptively selected covariate degree: %d\n',...
            %covariatedeg);
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
if isempty(targetlabels)
    coninds = 1:nlabels;
else
    coninds = find(strcmp(designvol.desc.features.unique.labels,...
        targetlabels));
end

if ~isempty(ignorelabels)
    ignoreinds = find(strcmp(designvol.desc.features.unique.labels,...
        ignorelabels));
    coninds = setdiff(coninds,ignoreinds,'stable');
end

assert(~isempty(coninds),...
    'found no labels matching targetlabels/ignorelabels %s/%s',...
    targetlabels,ignorelabels);
ncon = length(coninds);
%fprintf('storing estimates for %d conditions\n',ncon);

designvol = designvol(:,coninds);
predictornames = designvol.meta.features.labels;

% empty cells don't get read properly
if isempty(glmvarargs)
    glmvarargs = {};
end

% split the volumes
[designcell,epicell] = splitvol(split,designvol,epivol);
nsplit = numel(designcell);

glmcell = cell(nsplit,1);
nancell = cell(nsplit,1);
for s = 1:nsplit
    %fprintf('fitting split %d of %d...\n',s,nsplit);
    % split-specific data
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
    glmcell{s} = vol2glm(splitdesign,splitepi,glmclass,...
        covariatedeg,glmvarargs{:});
    nancell{s} = nanmask;
end

if usegpu
    for s = 1:nsplit
        glmcell{s} = structdata2class(glmcell{s},'gpuArray');
    end
end
