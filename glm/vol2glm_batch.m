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

% empty cells don't get read properly
if isempty(glmvarargs)
    glmvarargs = {};
end

% split the volumes
[designcell,epicell] = splitvol_batch(split,designvol,epivol,...
    'sgolayK',sgolayK,'sgolayF',sgolayF,'split',split,'targetlabels',...
    targetlabels,'ignorelabels',ignorelabels);
nsplit = numel(designcell);

% find a reasonable polynomial degree
if strcmp(covariatedeg,'adaptive')
    covariatedeg = vol2covdeg(epivol);
end

assert(isempty(sgolayK) || isempty(covariatedeg),['probably should not ' ...
    'detrend and include trend covariates at the same time']);

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
