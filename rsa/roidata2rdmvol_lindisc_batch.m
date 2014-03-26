% Linear discriminant-based RDM calculation. Given a set of ROIs
% (searchlight spheres or masks) and a pilab GLM instance, fit a linear
% discriminant for each ROI and compute some discriminant dissimilarity
% measure. 
%
% Historically, this used to wrap roidata2rdmvol_lindisc, but it
% now uses Process-based compute instead.
%
% INPUTS:
% rois: some MriVolume instance with rois in the sample dimension (rows)
% designvol: some BaseVolume instance (the RDM will have
%   designvol.nfeatures entries)
% epivol: some MriVolume instance (must be in register with rois)
%
% OUTPUTS:
% disvol: a volume of vectorised RDMs.
% splitdisvolcell: one disvol for each unique entry in split. Note that
%    compute is faster if you omit this output so only specify 2 nargout if
%    you really want to use this output.
%
% NAMED VARARGIN:
% glmclass: char defining GLM class (default GLM)
% glmvarargs: any additional arguments for GLM (e.g. k for RidgeGLM)
% crossvalidate: get split data discriminant (default true)
% cvsplit: indices to define GLM cvgroups (if crossvalidate==true)
% split: indices to define session-specific RDM compute. One RDM is
%   calculated for each unique entry in split and the resulting
%   session-specific RDMs are averaged.
% sterrunits: scale by standard error (default true)
% searchvol: if true, the resulting disvols are MriVolume. For this,
%   rois.nsamples must equal rois.nfeatures and epivol.nfeatures (default
%   false)
%
% [disvol,splitdisvolcell] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)
function [disvol,splitdisvolcell] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)

getArgs(varargin,{'split',[],'glmvarargs',{{}},'cvsplit',[],...
    'glmclass','GLM','sterrunits',1,'crossvalidate',1,...
    'minvoxeln',1,'searchvol',false,'crosscon',[],'minn',0});

% construct the RDM analysis
glman_rdm = rdms_lindisc_configureprocess('glmclass',glmclass,...
    'glmvarargs',...
    glmvarargs,'cvsplit',cvsplit,'sterrunits',sterrunits,...
    'crossvalidate',crossvalidate,'crosscon',crosscon,'ncon',...
    designvol.nfeatures,'setclass',class(epivol.data));

assert(isequal(epivol.xyz,rois.xyz),...
    'mismatched roi/epivol');
assert(isequal(epivol.meta.samples.chunks,...
    designvol.meta.samples.chunks),...
    'mismatched epivol/designvol');

if ~matlabpool('size')
    runfun = 'runrois_serial';
else
    runfun = 'runrois_spmd';
end


if isempty(split)
    split = ones(epivol.desc.samples.nunique.chunks,1);
end

if nargout>1
    % split the data into appropriately pre-processed cell arrays
    [designcell,epicell] = splitvol(split,designvol,epivol);
    nsplit = length(designcell);
    splitcell = cell(nsplit,1);
    % track NaN features - may appear in different runs if nan masking
    nanmask = false([nsplit rois.nsamples]);
    sl = ROIProcessor(rois,glman_rdm,minn,runfun);
    for sp = 1:nsplit
        fprintf('running rois for split %d of %d... ',sp,nsplit);
        tstart = clock;
        splitcell{sp} = call(sl,designcell{sp}.data,epicell{sp}.data,...
            epicoll{sp}.meta.samples.chunks);
        fprintf('finished in %.3fs\n',etime(clock,tstart));
        % track NaN features across splits may appear in different runs if
        % nan masking - nans should only really happen here if you enter
        % empty ROIs
        nanmask(sp,:) = any(isnan(splitcell{sp}),1);
    end % sp = 1:nsplit
    meanresult = matmean(splitcell{:});
    % remove any nan features from all sessdisvols
    anynan = any(nanmask,1);
    % make sure nans are consistent across sessions
    meanresult(:,anynan) = NaN;
else
    % can just build everything into one processor. This tends to be faster
    % because it puts more processing inside each parfor job.
    glman_rdm_split = SessionProcessor(split,glman_rdm);
    sl = ROIProcessor(rois,glman_rdm_split,minn,runfun);
    fprintf('running all rois... ')
    tstart = clock;
    meanresult = call(sl,designvol.data,epivol.data,...
        epivol.meta.samples.chunks);
    fprintf('finished in %.3fs\n',etime(clock,tstart));
end

% create output volume(s)
disvol = result2roivol(sl,meanresult,searchvol);

if nargout>1
    for sp = 1:nsplit
        % make sure nans are consistent
        splitdisvolcell{sp}(:,anynan) = NaN;
        splitdisvolcell{sp} = result2roivol(sl,splitcell{sp},searchvol);
    end
end
