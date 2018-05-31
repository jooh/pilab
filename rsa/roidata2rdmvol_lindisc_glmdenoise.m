% Linear discriminant-based RDM calculation with GLMDenoise fit. Given a
% set of ROIs (searchlight spheres or masks) and a pilab GLM instance, fit
% a linear discriminant for each ROI and compute some discriminant
% dissimilarity measure. 
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
% [disvol,denresult] = roidata2rdmvol_lindisc_glmdenoise(rois,designvol,epivol,varargin)
function [disvol,denresult] = roidata2rdmvol_lindisc_glmdenoise(rois,designvol,epivol,varargin)

getArgs(varargin,{'glmvarargs',{},'denoisetarget','raw',...
    'glmclass','GLM','sterrunits',0,...
    'minvoxeln',1,'searchvol',false,'crosscon',[],'minn',0,...
    'demean',false,'nperms',0,'covariatedeg','adaptive',...
    'numpcstotry',20,'pcR2cutoff',0,'pcstop',1.05,'pcR2cutoffmask',1,...
    'validationmode',false});

% for now we're going to fix this to strict leave-one-out for simplicity.
cvsplit = 1:epivol.desc.samples.nunique.chunks;

assert(nperms<2,'permutation testing is currently unsupported');

logstr('denoisetarget: %s\n',denoisetarget);

% set up the CV split loop - from searchselect
usplit = unique(cvsplit);
nsplit = length(usplit);

% get the covariates for each run
if strcmp(covariatedeg,'adaptive')
    covariatedeg = vol2covdeg(epivol);
end

glman_rdm = rdms_lindisc_configureprocess('glmclass',glmclass,...
    'glmvarargs',glmvarargs,'cvsplit',[1 2],'sterrunits',sterrunits,...
    'crossvalidate=1','crosscon',crosscon,'ncon',...
    designvol.nfeatures,'setclass',class(epivol.data),'demean',demean,...
    'nperms',nperms,'onewayvalidation',true);

% eventually we may have to be smarter about this - if we support
% permutation testing this will need to be switched by that too
% We cannot afford to do permutation testing with GLMdenoise. But the
% permutation test will be optimistic if we permute the data and the noise
% PCs. One solution is to use the CovariateGLM instance instead. This way
% the covariates would be left in place as we permute.
if ~matlabpool('size')
    runfun = 'runrois_serial';
else
    runfun = 'runrois_spmd';
end

chunks = epivol.meta.samples.chunks;
uchunks = epivol.desc.samples.unique.chunks;
nuchunks = epivol.desc.samples.nunique.chunks;
epimat = epivol.data;
designmat = designvol.data;

% prepare a detrended-only model to save a few cycles during
% crossvalidation later
for c = 1:nuchunks
    chunkind = chunks == uchunks(c);
    covariates{c} = cast(constructpolynomialmatrix(sum(chunkind),...
        0:covariatedeg),'like',epimat);
    detrendmodel(c) = CovariateGLM(designmat(chunkind,:),...
        epimat(chunkind,:),covariates{c});
end

splitres = cell(1,nsplit);
for sp = 1:nsplit
    logstr('running split %d of %d...\n',sp,nsplit);
    tstart = clock;
    thissplit = usplit(sp);
    trainind = cvsplit~=thissplit;

    % Noise pool processing to obtain noise covariates
    % figure out the noise pool (voxels with cv r2 < 0) for this training
    % split
    % there are various custom inputs we could consider but for now let's
    % keep it stock
    [noisepcs,denresult(sp)] = glm2denoisepcs(detrendmodel,trainind);

    % build a model for the training data (uchunks(trainind))
    switch denoisetarget
        case {'raw','mixandmatch'}
            % stock denoise - need DenoiseGLM instance
            denoisemodel = arrayfun(@(x)DenoiseGLM(detrendmodel(x).X,...
                detrendmodel(x).data,covariates{x},...
                noisepcs{x}),uchunks(trainind),'uniformoutput',0);
            denoisemodel = vertcat(denoisemodel{:});
            tuneprop = 'nnoisetouse';
            tunevalues = 0:numpcstotry;
        case 'denoised'
            % remove noise PCs from predicted data - just use noisepcs as
            % another covariate in CovariateGLM
            denoisemodel = arrayfun(@(x)CovariateGLM(detrendmodel(x).X,...
                detrendmodel(x).data,[covariates{x},noisepcs{x}]),...
                uchunks(trainind),'uniformoutput',0);
            denoisemodel = vertcat(denoisemodel{:});
            tuneprop = 'ncovtouse';
            % offset tune values since we don't want to tune the polynomial
            % degree...
            tunevalues = (0:numpcstotry) + covariatedeg+1;
        otherwise
            error('unknown denoisetarget: %s',denoisetarget);
    end
    % find the number of noise PCs in the training split
    [~,pcr2] = crossvalidateproperty(denoisemodel,tuneprop,tunevalues,...
        'rsquare');

    % use KK selection method - any voxel with r2>threshold in any analysis
    r2ok = any(pcr2>pcR2cutoff,1) & pcR2cutoffmask;
    % choose number of PCs
    chosen = NaN;
    % set 0 to no noise PC case.
    xvaltrend = median(pcr2(:,r2ok),2);
    curve = xvaltrend - xvaltrend(1);
    mx = max(curve);                   % store the maximum of the curve
    best = -Inf;                       % initialize (this will hold the best performance observed thus far)
    for p=0:numpcstotry
        % if better than best so far
        if curve(1+p) > best
            % record this number of PCs as the best
            chosen(sp) = p;
            best = curve(1+p);
            % if we are within pcstop of the max, then we stop.
            if best*pcstop >= mx
                break;
            end

        end
    end
    % you'll crash here e.g. if your data contains NaNs
    assert(~isnan(chosen(sp)),'unable to choose PC. model fit problem?');
    logstr('chosen number of noise PCs: %d\n',chosen(sp));

    % filter the dataset based on the chosen noise model
    [cvepi,cvrois] = intersectvols(epivol,rois);
    cvepimat = cvepi.data;
    cvdesignmat = designvol.data;
    cvchunk = zeros(epivol.nsamples,1);
    % denoise the training data/design
    for c = asrow(uchunks(trainind))
        chunkind = chunks==c;
        pmat = projectionmatrix([covariates{c} ...
            noisepcs{c}(:,1:chosen(sp))]);
        cvepimat(chunkind,:) = pmat * cvepimat(chunkind,:);
        cvdesignmat(chunkind,:) = pmat * cvdesignmat(chunkind,:);
        % all training runs get put in the same chunk
        cvchunk(chunkind) = 1;
    end
    
    % handle the test chunk
    chunkind = chunks==thissplit;
    assert(all(cvchunk(chunkind)==0),'non-independent split detected');
    cvchunk(chunkind) = 2;
    switch denoisetarget
        case 'raw'
            pmat = projectionmatrix(covariates{thissplit});
        case {'denoised','mixandmatch'}
            pmat = projectionmatrix([covariates{thissplit} ...
                noisepcs{thissplit}(:,1:chosen(sp))]);
        otherwise
            error('unknown denoisetarget: %s',denoisetarget);
    end
    cvepimat(chunkind,:) = pmat * cvepimat(chunkind,:);
    cvdesignmat(chunkind,:) = pmat * cvdesignmat(chunkind,:);

    roiiterator = ROIProcessor(cvrois,glman_rdm,minn,runfun);
    splitres{sp} = call(roiiterator,cvepimat,cvdesignmat,cvchunk);
    logstr('finished in %s\n',seconds2str(etime(clock,tstart)));

    if validationmode
        tdnstart = clock;
        designcell = {denoisemodel.X};
        datacell = cellfun(@transpose,{denoisemodel.data},...
            'uniformoutput',0);
        result = GLMdenoisedata(designcell,datacell,NaN,epivol.frameperiod,...
            'assumeconvolved',NaN,struct('numboots',0,'maxpolydeg',4),[]);
        logstr('finished GLMdenoise in %s\n',seconds2str(etime(clock,...
            tdnstart)));
    end
end
% 4) average the splits.
meanresult = matmean(splitres{:});

% combine the results
chosen = num2cell(chosen);
[denresult.chosen] = deal(chosen{:});

% create output volume(s)
disvol = result2roivol(roiiterator,meanresult,searchvol);
