%
%
% (varargin,{'split',[],'glmvarargs',{{}},'cvsplit',[],...
    %'glmclass','GLM','sterrunits',1,'crossvalidate',1,...
    %'minvoxeln',1,'crosscon',[],...
    %'subsplit',[],'maskns',[],'predictor',[],...
    %'rsaclass','RankRSA','rsaclassargs',{},'minn',0});
%
% [meandisvol,sessdisvolcell,roispheres] = roidata2rdmvol_lindisc_searchselect(rois,designvol,epivol,masks,varargin)
function [meandisvol,sessdisvolcell,roispheres] = roidata2rdmvol_lindisc_searchselect(rois,designvol,epivol,masks,varargin)

getArgs(varargin,{'split',[],'glmvarargs',{},'cvsplit',[],...
    'glmclass','GLM','sterrunits',1,'crossvalidate',1,...
    'minvoxeln',1,'crosscon',[],...
    'subsplit',[],'maskns',[],'predictor',[],...
    'rsaclass','RankRSA','rsaclassargs',{},'minn',0});

% configure the basic RDM analysis
glman_rdm = rdms_lindisc_configureprocess('glmclass',glmclass,...
    'glmvarargs',...
    glmvarargs,'cvsplit',cvsplit,'sterrunits',sterrunits,...
    'crossvalidate',crossvalidate,'crosscon',crosscon,'ncon',...
    designvol.nfeatures,'setclass',class(epivol.data));

if ~matlabpool('size')
    runfun = 'runrois_serial';
else
    runfun = 'runrois_spmd';
end

if ~iscell(rsaclassargs)
    rsaclassargs = {rsaclassargs};
end

% configure rsa
rsaconstruct = RSAConstructor(rsaclass,predictor,rsaclassargs{:});
rsafit = GLMProcessor('fit',[],1);
rsafull = RSAMetaProcessor(rsaconstruct,rsafit);
fullprocessor = SecondLevelProcessor(glman_rdm,rsafull);

nsizes = numel(maskns);

% RSA contrast matrix
% TODO how deal with cross-decoding here? Probably need to set up
% processors here. 
% the number of processors in glman_rdm gives the number of inputs to
% provide to discriminant (and the number of outputs to expect). 
% this is a little hacky but for now it's a workable way to retrieve the
% training conmats
nweight = numel(glman_rdm.processor);
conmats = arrayfun(@(thisproc)thisproc.varargs{4},glman_rdm.processor,...
    'uniformoutput',false);
discriminator = GLMProcessor('discriminant',[],nweight,conmats{:});
% pass handle to same constructor as main glman_rdm
getweights = GLMMetaProcessor(glman_rdm.constructor,discriminator);
npairs = size(conmats{1},1);

% there are 3 split levels here: 
% split: defines which sessions are used for computing RDMs (the
% final RDM for the searchlight / ROI is the average over all
% splits)
% subsplit: defines which part of each split is used for
% searchlight mapping and which is used to generate the final
% independent RDM
% cvsplit: defines how the LDt RDM is cross-validated within each
% session
if ischar(subsplit)
    subsplit = eval(subsplit);
end
usub = unique(subsplit);
nsub = length(usub);

if isempty(split)
    split = ones(designvol.desc.samples.nunique.chunks,1);
end
nsplit = numel(unique(split));
if nsplit > 1
    [designcell,epicell] = splitvol(split,designvol,epivol);
else
    designcell = {designvol};
    epicell = {epivol};
end

if ~iscell(masks)
    masks = repmat({masks},[1 nsplit]);
else
    assert(numel(masks)==nsplit,...
        'number of masks must match number of splits');
end

% different logic here now - do the entire processing stream for each
% split.
for sp = 1:nsplit
    fprintf('running searchselect for split %d of %d...\n',sp,nsplit);
    tstart = clock;
    epivol = epicell{sp};
    designvol = designcell{sp};
    % different mask for each subsp split, potentially
    maskvol = masks{sp};
    % make sure we have the same masks, ROIs and voxels across splits
    [rois,epivol,maskvol] = intersectvols(rois,epivol,maskvol);
    uchunks = epivol.desc.samples.unique.chunks;

    % for each subsp-split (roi defining and testing split)
    for subsp = 1:nsub
        fprintf('running sub-split %d of %d...\n',subsp,nsub)
        % figure out which data will be used to identify ROI and
        % which to test
        thissuper = usub(subsp);

        % for each mask
        for mask = 1:maskvol.nsamples
            maskname = maskvol.meta.samples.names{mask};
            maskind = maskvol.data(mask,:)~=0;
            epimask = epivol(:,maskind);
            roimask = rois(maskind,maskind);

            % compute searchlight RDMs in training split (separately
            % for different sessions if desired)
            trainind = subsplit~=thissuper;
            trainchunks = uchunks(trainind);
            epitrain = epimask.selectbymeta('chunks',trainchunks);
            designtrain = designvol.selectbymeta('chunks',trainchunks);
            sl = ROIProcessor(roimask,fullprocessor,minn,runfun);
            % produces equivalent output to using
            % roidata2rdmvol_lindisc_batch and plugging the resulting RDMs
            % into roidata_rsa
            r = sl.call(designtrain.data,epitrain.data,...
                epitrain.meta.samples.chunks);

            % make a new roivol containing ROIs generating according each
            % roi size
            roimat = false([nsizes epimask.nfeatures]);
            definingcon = repmat({predictor.name},[nsizes 1]);
            names = cell(nsizes,1);
            thresholds = NaN([nsizes,1]);
            weights = cell(nsizes,nweight);
            for n = 1:nsizes
                thissize = maskns(n);
                % make the ROI volume in the test data according to
                % the ranked effect in train data
                [roimat(n,:),thresholds(n)] = ...
                    selectbysearchlight(roimask,r,thissize,'selectmode',...
                    'union');
                names{n} = sprintf('%s %s (%d voxels)',...
                    maskname,predictor.name,thissize);
                % obtain some set of weights
                [weights{n,1:nweight}] = call(getweights,...
                    designtrain.data,epitrain.data(:,roimat(n,:)~=0),...
                    designtrain.meta.samples.chunks);
            end % n = 1:nsizes

            roispheres.(maskname){sp,subsp} = MriVolume(roimat,epimask,...
                'metasamples',struct('names',{names},...
                'thresholds',thresholds,'definingcon',{definingcon},...
                'nfeatures_target',maskns,'masknames',...
                {repmat({maskname},[nsizes 1])}));

            % done with training data - clear here to be extra safe
            % with independence
            clear epitrain designtrain trainind trainchunks trainmodel sl

            % TEST DATA PROCESSING
            testind = subsplit==thissuper;
            testchunks = uchunks(testind);
            epitest = epimask.selectbymeta('chunks',testchunks);
            designtest = designvol.selectbymeta('chunks',testchunks);
            % now we just need to apply the test method with the weights.
            % this involves iteration so I suppose another ROI Processor is
            % in order technically (but we'd struggle to pass the weights
            % appropriately so let's not)
            testinfot = NaN(npairs,nsizes,nweight);
            for n = 1:nsizes
                testmodel = GLM(designtest.data,...
                    epitest.data(:,roimat(n,:)~=0));
                % now the fun twist is that there may be multiple weights
                for w = 1:nweight
                    testinfot(:,n,w) = infotmap(testmodel,weights{n,w},...
                        conmats{w});
                end
            end
            % collapse weight dimension to average cross-decoders (e.g., A
            % to B and B to A)
            testinfot = mean(testinfot,3);
            roip = ROIProcessor(roispheres.(maskname){sp,subsp},[]);
            splitdisvol{mask,subsp} = result2roivol(roip,testinfot);

            % clear for independence safety
            clear epitest designtest testind testchunks testmodel roip
        end % mask = 1:maskvol.nsamples 

    end % subsp = 1:nsub

    % combine the sub-splits into one session disvol
    for mask = 1:maskvol.nsamples
        for subsp = 1:nsub
            sessdisvolcell{sp} = cat(2,splitdisvol{:,subsp});
        end
    end
    fprintf('finished split in %s\n',seconds2str(etime(clock,tstart)));
end % sp = 1:nsplit

if nsplit<2
    % quick and easy
    meandisvol = sessdisvolcell{1};
else
    % average into one meandisvol
    for sp = 1:nsplit
        disdata{sp} = sessdisvolcell{sp}.data;
        disfeat(sp) = sessdisvolcell{sp}.meta.features;
    end
    meandis = matmean(disdata{:});
    % combine the meta data from all splits for roi and disvol
    roimeta = collapsestruct(disfeat,@mean);
    meandisvol = BaseVolume(meandis,...
        'metafeatures',roimeta);
end
