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
% batchsize: number of ROIs to run in one call to parfor (default 1000)
% searchvol: if true, the resulting disvols are MriVolume. For this,
%   rois.nsamples must equal rois.nfeatures and epivol.nfeatures (default
%   false)
%
% [disvol,splitdisvolcell] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)
function [disvol,splitdisvolcell] = roidata2rdmvol_lindisc_batch(rois,designvol,epivol,varargin)

% TODO - cross-class RDM support how? Need some kind of index specifying
% a and b.
getArgs(varargin,{'split',[],'glmvarargs',{{}},'cvsplit',[],...
    'glmclass','GLM','sterrunits',1,'crossvalidate',1,...
    'minvoxeln',1,'batchsize',10000,'searchvol',false,'crossinds',[]});

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

if ischar(crossinds)
    crossinds = eval(crossinds);
end

if sterrunits
    testmeth = 'infotmap';
else
    testmeth = 'infomahalanobis';
end

% assemble processors
glmspec = GLMConstructor('GLM',cvsplit);

if isempty(crossinds)
    % straight RDM
    np = nchoosek(designvol.nfeatures,2);
    % contrast vector of correct class
    cons = allpairwisecontrasts(feval(class(epivol.data),...
        designvol.nfeatures));
    if crossvalidate
        rdm = GLMProcessor('cvclassificationrun',[],1,'discriminant',...
            testmeth,[np 1],cons);
        % only necessary if CV
        glman_rdm = GLMMetaProcessor(glmspec,rdm,@(x)mean(x,3));
    else
        error('currently unsupported');
    end
else
    % cross-class RDM - first set up contrast vectors
    assert(crossvalidate,'must crossvalidate if crossinds are present');
    nc = numel(crossinds{1});
    assert(nc==numel(crossinds{2}),'mismatched crossinds');
    assert(nc<=(designvol.nfeatures/2),...
        'bad crossinds for designvol.nfeatures');
    np = nchoosek(nc,2);
    cons = allpairwisecontrasts(feval(class(epivol.data),nc));
    fullmat = zeros(np,designvol.nfeatures,class(cons));
    crosscon{1} = fullmat;
    crosscon{1}(:,crossinds{1}) = cons;
    crosscon{2} = fullmat;
    crosscon{2}(:,crossinds{2}) = cons;
    % then processors
    rdmcross(1) = GLMProcessor('cvcrossclassificationrun',[],1,...
        'discriminant',testmeth,[np 1],crosscon{1},crosscon{2});
    rdmcross(2) = GLMProcessor('cvcrossclassificationrun',[],1,...
        'discriminant',testmeth,[np 1],crosscon{2},crosscon{1});
    glman_rdm = GLMMetaProcessor(glmspec,rdmcross,@(x)mean(x,3));
end


if nargout>1
    % split the data into appropriately pre-processed cell arrays
    [designcell,epicell] = splitvol(split,designvol,epivol);
    nsplit = length(designcell);
    splitcell = cell(nsplit,1);
    % track NaN features - may appear in different runs if nan masking
    nanmask = false([nsplit rois.nsamples]);
    sl = ROIProcessor(rois,glman_rdm,batchsize);
    for sp = 1:nsplit
        fprintf('running rois for split %d of %d...\n',sp,nsplit);
        tstart = clock;
        splitcell{sp} = call(sl,designcell{sp},epicell{sp});
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
    sl = ROIProcessor(rois,glman_rdm_split,batchsize);
    meanresult = sl.call(designvol,epivol);
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
