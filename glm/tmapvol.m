% Compute T maps against baseline for an epi volume and a design volume
% tvol = tmapvol(designvol,epivol,[inds])
%
% inds are indices into uniquelabels (defaults to 1:nlabels)
function tvol = tmapvol(designvol,epivol,inds,covariatedeg)

if nargin<3
    inds = 1:designvol.desc.features.nunique.labels;
end
nt = length(inds);

if strcmp(covariatedeg,'adaptive')
    % Kay's rule for deciding on n covariates based on run duration.
    covariatedeg = round(epivol.nsamples*epivol.frameperiod/60/2);
    assert(~isempty(covariatedeg),'failed to identify covariatedeg');
end

% fit
fit = CovGLM(designvol.data,epivol.data,covariatedeg);
% prepare output matrix
tdata = NaN([nt epivol.nfeatures]);
% concessions to parfor to avoid memory overhead (actually, further testing
% suggests you get better performance without parfor here, so disabled for
% now)
colinds = 1:epivol.nfeatures;
labelinds = designvol.desc.features.inds.labels;
regs = zeros(1,fit.npredictors);
% iterate to make t maps
for t = 1:nt
    % make the vector for this label (considering any extra nuisance vars)
    cv = regs;
    cv(labelinds==inds(t)) = 1;
    % mean over standard error
    tdata(t,colinds) = fit.tmap(cv);
end

% make the Volume instance
tvol = MriVolume(tdata,epivol,'metasamples',struct('labels',...
    {designvol.desc.features.unique.labels(inds)'}));
