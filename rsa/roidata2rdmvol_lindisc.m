% given a set of ROIs (searchlight spheres or masks) and a pilab GLM
% instance, fit a linear discriminant for each ROI and compute some
% discriminant dissimilarity measure according to parameters in the ts
% struct.
%
% roidata2rdmvol is used when you have one response pattern per condition
% and want some pdist-based distance metric on the patterns. This function
% is used when you have raw data and a design matrix and want to take
% advantage of the better error covariance matrix estimated from the full
% model residual or to cross-validate a split-data RDM, e.g. for LDAt RDMs.
%
% named varargs:
%
% glmclass: char defining GLM class (default GLM)
% glmvarargs: any additional arguments for GLM (e.g. k for RidgeGLM)
% split: indices to define GLM cvgroups (if crossvalidate) - NB NOT used in
%   vol2glm_batch
% sterrunits: scale by standard error (default false)
% crossvalidate: get split data discriminant (default false)
% batchsize: number of ROIs to run in one call to parfor (default 5000)
% 
% The gist of these options are:
% ~split && ~sterrunits = mahalanobis distance
% ~split && sterrunits = 'whitened euclidean' (?)
% split && ~sterrunits = mahalanobis classifier distance
% split && sterrunits = LDAt RDM
%
% disvol = roidata2rdmvol_lindisc(rois,designvol,epivol,[varargin])
function disvol = roidata2rdmvol_lindisc(rois,designvol,epivol,varargin)

ts = varargs2structfields(varargin,struct(...
    'split',[],...
    'glmclass','GLM','glmvarargs',{},'sterrunits',false,'crossvalidate',...
    false,'minvoxeln',0,'batchsize',5000));

if ~iscell(ts.split)
    % so we can easily deal to cvgroup field later
    ts.split = num2cell(ts.split);
end

% preallocate result (dissimilarity by roi)
dissimilarities = NaN([nchoosek(designvol.nfeatures,2) rois.nsamples],...
    class(epivol.data));

% check for nans
nanmask = ~any(isnan(epivol.data),1);

% pairwise contrasts
conmat = feval(class(epivol.data),allpairwisecontrasts(designvol.nfeatures));

if ts.sterrunits
    testmeth = 'infotmap';
else
    testmeth = 'infomahalanobis';
end

% pad out last batch with NaNs
batchinds = [1:rois.nsamples ...
    NaN([1 ts.batchsize-rem(rois.nsamples,ts.batchsize)])];
batchmat = reshape(batchinds,...
    [ts.batchsize ceil(rois.nsamples/ts.batchsize)]);
nbatch = size(batchmat,2);

assert(isequal(epivol.xyz,rois.xyz),'mismatched roi/epivol');
assert(isequal(epivol.meta.samples.chunks,designvol.meta.samples.chunks),...
    'mismatched epivol/designvol');
% now store just array data for extra speed
designmat = designvol.data;

% batch out ROIs to allow a smaller epivol. Avoids Matlab memory problems
% when parfor involves > 2GB of data and avoids passing a huge epivol
% around when only a small part of it will actually be used.
for batch = 1:nbatch
    fprintf('running RDMs for batch %d of %d...\n',batch,nbatch);
    tic;
    if any(isnan(batchmat(:,batch)))
        thisbatchsize = find(isnan(batchmat(:,batch)),1,'first')-1;
    else
        thisbatchsize = ts.batchsize;
    end
    % voxels in any roi and not nan
    batchvox = any(full(rois.data(batchmat(1:thisbatchsize,batch),:)~=0)...
        ,1) & nanmask;
    % pick these for batch-specific epivol
    % may need linind2featind here
    batchepis = epivol(:,batchvox);
    batchrois = rois(:,batchvox);
    assert(numel(batchvox)==epivol.nfeatures,'mismatched mask and data');
    assert(numel(batchvox)==rois.nfeatures,'mismatched mask and data');

    % just arrays for speed
    chunks = epivol.meta.samples.chunks;
    roimat = batchrois.data;
    epimat = batchepis.data;

    % now we should just be able to run parfor directly
    if ts.crossvalidate
        parfor b = batchmat(1:thisbatchsize,batch)'
            % skip empty rois (these come out as NaN)
            validvox = full(roimat(batchinds(b),:)~=0);
            if ~any(validvox) || sum(validvox)<ts.minvoxeln
                % empty or too small roi
                continue
            end
            thismodel = array2glm(designmat,epimat(:,validvox),chunks,...
                ts.glmclass,ts.glmvarargs{:});
            % split defines crossvalidation split in GLM (NB in other contexts
            % split may get passed to vol2glm_batch instead to make one GLM
            % instance per split).
            [thismodel.cvgroup] = ts.split{:};
            cvres = cvclassificationrun(thismodel,'discriminant',testmeth,...
                [],conmat);
            % result - mean across splits
            dissimilarities(:,b) = mean(cvres,3);
        end
    else
        parfor b = batchmat(1:thisbatchsize,batch)'
            % skip empty rois (these come out as NaN)
            validvox = full(roimat(batchinds(b),:)~=0);
            if ~any(validvox) || sum(validvox)<ts.minvoxeln
                % empty or too small roi
                continue
            end
            thismodel = array2glm(designmat,epimat(:,validvox),chunks,...
                ts.glmclass,ts.glmvarargs{:});
            % just self-fit
            w = discriminant(thismodel,conmat);
            dissimilarities(:,b) = feval(testmeth,thismodel,...
                w,conmat);
        end
    end
    fprintf('finished batch in %s\n',seconds2str(toc));
end

% convert to volume - here it is a problem that the result may have
% different nfeatures than the mask (e.g. for ROI analysis or when we do
% not run all possible searchlight spheres)
if rois.nsamples == rois.nfeatures
    % simple case - assume that samples and features are in register
    disvol = MriVolume(dissimilarities,rois,'metafeatures',struct(...
        'names',{rois.meta.samples.names}));
else
    % complicated case - need to forget the mask and write out a mask-less
    % volume. But save coordinates of ROIs to enable sanity checks later
    coords = cell(1,rois.nsamples);
    nvox = NaN([1 rois.nsamples]);
    for c = 1:rois.nsamples
        % compute centre of mass for this ROI
        coords{c} = round(mean(rois.linind2coord(rois.linind(...
            rois.data(c,:)~=0)),2));
        nvox(c) = sum(rois.data(c,:)~=0);
    end
    % make a mask-less volume 
    disvol = MriVolume(dissimilarities,[],'metafeatures',struct(...
        'names',{rois.meta.samples.names'},'centreofmass',{coords},...
        'nfeatures',nvox),'header',rois.header);
end
