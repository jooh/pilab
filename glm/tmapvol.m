% Compute T maps against baseline for an epi volume and a design volume
% tvol = tmapvol(designvol,epivol,[inds])
%
% inds are indices into uniquelabels (defaults to 1:nlabels)
function tvol = tmapvol(designvol,epivol,inds)

if nargin<3
    inds = 1:designvol.desc.features.nunique.labels;
end
nt = length(inds);

% fit
model = GLM(designvol.data,epivol.data);
% prepare output matrix
tdata = NaN([nt epivol.nfeatures]);
% concessions to parfor to avoid memory overhead (actually, further testing
% suggests you get better performance without parfor here, so disabled for
% now)
colinds = 1:epivol.nfeatures;
labelinds = designvol.desc.features.inds.labels;
% iterate to make t maps
for t = 1:nt
    % make the vector for this label
    cv = double(labelinds==inds(t));
    % mean over standard error
    tdata(t,colinds) = model.tmap(cv);
end

% make the Volume instance
tvol = MriVolume(tdata,epivol,'metasamples',struct('labels',...
    {designvol.desc.features.unique.labels(inds)'}));
