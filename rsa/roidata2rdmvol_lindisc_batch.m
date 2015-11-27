% Linear discriminant-based RDM calculation. Given a set of ROIs
% (searchlight spheres or masks) and a pilab GLM instance, fit a linear
% discriminant for each ROI and compute some discriminant dissimilarity
% measure. 
%
% Historically, this used to wrap roidata2rdmvol_lindisc, but it
% now uses Process-based compute instead.
%
% INPUTS:
% rois: some SPMVolume instance with rois in the sample dimension (rows)
% designvol: some Volume instance (the RDM will have
%   designvol.nfeatures entries)
% epivol: some SPMVolume instance (must be in register with rois)
%
% OUTPUTS:
% disvol: a volume of vectorised RDMs.
% splitdisvolcell: one disvol for each unique entry in split. Note that
%    compute is faster if you omit this output so only specify 2 nargout if
%    you really want to use this output.
% permres: a matrix of permuted results stacked in third dim (true perm is
%   first entry)
% sesspermres: cell array with permuted results for each session split.
%
% NAMED VARARGIN:
% glmclass: char defining GLM class (default GLM)
% glmvarargs: any additional arguments for GLM (e.g. k for RidgeGLM)
% crossvalidate: get split data discriminant (default true)
% cvsplit: indices to define GLM cvgroups (if crossvalidate==true)
% split: indices to define session-specific RDM compute. One RDM is
%   calculated for each unique entry in split and the resulting
%   session-specific RDMs are averaged.
% sterrunits: scale by standard error (default false)
% searchvol: if true, the resulting disvols are SPMVolume. For this,
%   rois.nsamples must equal rois.nfeatures and epivol.nfeatures (default
%   false)
% nperms: number of permutations to use for permuteruns (default 0)
%
% [disvol,splitdisvolcell,permres,sesspermres] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)
function [disvol,splitdisvolcell,permres,sesspermres] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)

getArgs(varargin,{'split',[],'glmvarargs',{},'cvsplit',[],...
    'glmclass','GLM','sterrunits',0,'crossvalidate',1,...
    'searchvol',false,'crosscon',[],'minn',0,...
    'demean',false,'nperms',0,'onewayvalidation',false});

% construct the RDM analysis
if isempty(cvsplit)
    cvsplit = 1:epivol.desc.samples.nunique.chunks;
    assert(isempty(split),'cvsplit must be provided if split is defined');
end

glman_rdm = rdms_lindisc_configureprocess('glmclass',glmclass,...
    'glmvarargs',...
    glmvarargs,'cvsplit',cvsplit,'sterrunits',sterrunits,...
    'crossvalidate',crossvalidate,'crosscon',crosscon,'ncon',...
    designvol.nfeatures,'setclass',class(epivol.data),'demean',demean,...
    'nperms',nperms,'onewayvalidation',onewayvalidation);

if isfield(epivol,'xyz')
    % support non-MriVolumes
    assert(isequal(epivol.xyz,rois.xyz),...
        'mismatched roi/epivol');
end
assert(isequal(epivol.meta.samples.chunks,...
    designvol.meta.samples.chunks),...
    'mismatched epivol/designvol');

% serial ROI iteration if no matlabpool OR if we are running a permutation
% test (since this already uses parfor).
% you'd think this would make very little difference the speedup is about
% 3.5x compared to running permutation tests with this setting still in
% runrois_spmd and letting Matlab sort out the parallelisation.
if ~matlabpool('size') || nperms > 1
    runfun = 'runrois_serial';
else
    runfun = 'runrois_spmd';
end
logstr('running ROIProcessing with %s\n',runfun);

if isempty(split)
    split = ones(epivol.desc.samples.nunique.chunks,1);
end

if nargout>1
    % split the data into appropriately pre-processed cell arrays
    [designcell,epicell] = splitvol(split,designvol,epivol);
    nsplit = length(designcell);
    sesspermres = cell(nsplit,1);
    % track NaN features - may appear in different runs if nan masking
    nanmask = cell(nsplit,1);
    sl = ROIProcessor(rois,glman_rdm,minn,runfun);
    for sp = 1:nsplit
        logstr('running rois for split %d of %d...\n',sp,nsplit);
        tstart = clock;
        sesspermres{sp} = call(sl,epicell{sp}.data,designcell{sp}.data,...
            epicell{sp}.meta.samples.chunks);
        logstr('finished in %s\n',seconds2str(etime(clock,tstart)));
        % track NaN features across splits may appear in different runs if
        % nan masking - nans should only really happen here if you enter
        % empty ROIs
        nanmask{sp} = any(isnan(sesspermres{sp}),1);
    end % sp = 1:nsplit
    permres = matmean(sesspermres{:});
    nanmask = cat(1,nanmask{:});
    % remove any nan features from all sessdisvols
    anynan = any(nanmask,3);
    anynan = any(anynan,1);
    % make sure nans are consistent across sessions
    permres(:,anynan,:) = NaN;
else
    % can just build everything into one processor. This tends to be faster
    % because it puts more processing inside each parfor job.
    glman_rdm_split = SessionProcessor(split,glman_rdm);
    sl = ROIProcessor(rois,glman_rdm_split,minn,runfun);
    logstr('running all rois... ')
    tstart = clock;
    permres = call(sl,epivol.data,designvol.data,...
        epivol.meta.samples.chunks);
    logstr('finished in %s\n',seconds2str(etime(clock,tstart)));
end

% deal with permutation test
meanresult = permres(:,:,1);

% create output volume(s)
disvol = result2roivol(sl,meanresult,searchvol);

if nargout>1
    for sp = 1:nsplit
        % make sure nans are consistent
        sesspermres{sp}(:,anynan,:) = NaN;
        splitdisvolcell{sp} = result2roivol(sl,sesspermres{sp}(:,:,1),...
            searchvol);
    end
end
