% batch convert designvol / epivol to a cell of entries instances. Used
% internally in vol2glm_batch
%
% named varargs:
%
% sgolayK: degree of optional Savitsky-Golay detrend 
% sgolayF: filter size of optional Savitsky-Golay detrend
% split: indices to define which runs go into which glmcell entry
% targetlabels: labels to explicitly include (default all)
% ignorelabels: labels to explicitly exclude (default none)
%
% [designcell,epicell] = splitvol_batch(split,designvol,epivol,varargin)
function [designcell,epicell] = splitvol_batch(split,designvol,epivol,varargin)

getArgs(varargin,{'sgolayK',[],'sgolayF',[],'targetlabels',{},...
    'ignorelabels',{}});

nchunks = epivol.desc.samples.nunique.chunks;

% detrend before splitting
if ~isempty(sgolayK)
    fprintf('detrending (K=%d,F=%d)\n',sgolayK,sgolayF);
    epivol.sgdetrend(sgolayK,sgolayF);
    designvol.sgdetrend(sgolayK,sgolayF);
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

if ~isequal(coninds,1:nlabels)
    designvol = designvol(:,coninds);
end
predictornames = designvol.meta.features.labels;

% split the volumes
[designcell,epicell] = splitvol(split,designvol,epivol);
