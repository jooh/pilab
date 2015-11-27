% do a random effects analysis on the data in struct array subres (one
% entry per subject, or a single entry 'groupres' where subjects are
% stacked in third dimension). Makes various assumptions about how subres
% is organised - see roidata_rsa and roidata_glm for reference.
%
% Named inputs:
% nperm: 1 
% nboot: 0
% targetfield: 'r'
% transfun: string that gets feval'ed or a function handle to apply some
%   transform to the data before doing the rfx analysis (e.g. 'atanh' for
%   fisher transform, or something like @(x)x-.5 to mean correct proportion
%   correct classification performance).
% assumeregister: false (if true, we just stack the single subject results,
% consequences be damned)
% varsmoothmask: []
% varsmoothfwhm: 8 (mm)
% keepnans: (false) preserve nan conditions in outputs
% minn: 2 skip any roi with fewer than this n
% multpdim: ('roi') can also be 'con' for FWE/FDR correction over conditions
%   instead
% contrasts: calculate differences between particular conditions. char that
%   gets evaled or struct with fields:
%       name
%       conplus
%       conminus
%       tail
% customfits: obtain a result by applying some function to the data. char
%   that gets fevaled or struct with the following mandatory fields:
%       name
%       funhand
%       cons
%       tail
%   the syntax for calling funhand is as follows, so the customfits can
%   contain custom settings that determine funhand's behavior:
%       groupres.(targetfield)(end+1,:,:) = customfits(c).funhand(...
%       groupres.(targetfield)(conind,:,:),customfits(c));
% customfun: catch-all custom processing. Gets called as 
%       meanres.custom{nc,r} = feval(customfun{nc},roidata(:,:,r),...
%           meanres.rows_contrast);
%
% [meanres,groupres,nulldist,bootdist] = roidata_rfx(subres,[varargin])
function [meanres,groupres,nulldist,bootdist] = roidata_rfx(subres,varargin)
    
getArgs(varargin,{'nperm',1,'nboot',0,'targetfield','r','transfun',[],...
    'assumeregister',false,'varsmoothmask',[],'varsmoothfwhm',8,...
    'contrasts',emptystruct('name','conplus','conminus','tail'),...
    'customfits',emptystruct('name','funhand','cons','tail'),'minn',2,...
    'keepnans',false,'multpdim','roi','customfun',[]});

% we assume that contrasts are identical.  we could support
% this case in theory but let's not for now since it likely
% indicates a user error.
concell = {subres.rows_contrast};
assert(isequal(concell{1},concell{:}),...
    'different predictors in different subjects');
ncon = length(concell{1});

if ~iscell(customfun)
    if isempty(customfun)
        customfun = {};
    else
        customfun = {customfun};
    end
end
ncustom_rfx = numel(customfun);

isgroupres = false;
nsub = length(subres);
targets = intersect(fieldnames(subres),...
    {targetfield,'pperm','ppara','medianboot','medianste',...
    'nfeatures','tail'});
if numel(subres)==1 && isfield(subres,'z_subject')
    % groupres mode
    isgroupres = true;
    groupres = subres;
    names = groupres.z_subject;
else
    names = {subres.name};
    if assumeregister
        % just take rois from first subject
        urois = subres(1).cols_roi;
    else
        % find all possible ROIs (not all subjects will have all ROIs
        % necessarily)
        allrois = {subres.cols_roi};
        allrois = horzcat(allrois{:});
        urois = unique(allrois);
    end
    nroi = length(urois);
    dat = NaN([ncon nroi nsub]);
    groupres = struct('rows_contrast',{subres(1).rows_contrast},...
        'cols_roi',{urois},targetfield,dat,'nfeatures',NaN([1 nroi nsub]),...
        'tail',{subres(1).tail});

    % populate the groupres struct
    if assumeregister
        % hey, this is easy
        groupres = collapsestruct(subres,@zcat);
        assert(~any(isnan(groupres.(targetfield)(:))),...
            'nans in targetfield not supported in assumeregister mode');
    else
        % have to be careful not to average over different rois
        for s = 1:nsub
            for r = 1:length(subres(s).cols_roi)
                thisroi = subres(s).cols_roi{r};
                indroi = strcmp(groupres.cols_roi,thisroi);
                for t = targets(:)'
                    tstr = t{1};
                    if strcmp(tstr,'tail')
                        assert(isequal(groupres.tail,subres.tail),...
                            'mismatched tails across subjects');
                    else
                        groupres.(tstr)(:,indroi,s) = ...
                            subres(s).(tstr)(:,r);
                    end
                end
            end
        end
    end
end

if isfield(subres,'custom')
    ncustom = numel(subres(1).custom);
    groupres.custom = cell(size(subres(1).custom));
    for s = 1:nsub
        assert(isequal(ncustom,numel(subres(s).custom)),...
            'mismatched numel(custom) across subjects');
        % work out the ROI indices for this subject
        [~,ginds,sinds] = intersect(groupres.cols_roi,subres(s).cols_roi,...
            'stable');
        for c = 1:ncustom
            assert(ismatrix(subres(s).custom{c}),['>2d custom fields '...
                'are not supported']);
            groupres.custom{c} = zcat(groupres.custom{c},...
                NaN([size(subres(s).custom{c},1),nroi]));
            % so we make no assumptions about the shape of custom here except
            % that ROIs are in second dim
            groupres.custom{c}(:,ginds,end) = subres(s).custom{c}(:,sinds);
        end
    end
end

if ~isempty(transfun)
    if ischar(transfun)
        transfun = str2func(transfun);
    end
    % apply some transform, e.g. fisher
    groupres.(targetfield) = transfun(groupres.(targetfield));
end
assert(isreal(groupres.(targetfield)),'imaginary output!');

tail = groupres.tail;

% add in custom fits
if ~isempty(customfits)
    if ischar(customfits)
        customfits = feval(customfits,groupres.rows_contrast);
    end
    assert(isstruct(customfits),'customfits must be a struct');
    for c = 1:length(customfits)
        if ischar(customfits(c).funhand)
            customfits(c).funhand = str2func(customfits(c).funhand);
        end
        conind = findconind(groupres.rows_contrast,customfits(c).cons);
        groupres.(targetfield)(end+1,:,:) = customfits(c).funhand(...
            groupres.(targetfield)(conind,:,:),customfits(c));
        groupres.rows_contrast{end+1} = customfits(c).name;
        tail{end+1} = customfits(c).tail;
    end
end

% add in the contrasts
if ~isempty(contrasts)
    if ischar(contrasts)
        contrasts = feval(contrasts,groupres.rows_contrast);
    end
    assert(isstruct(contrasts),'contrasts must be a struct');
    for c = 1:length(contrasts)
        % find the rows corresponding to each condition
        posind = findconind(groupres.rows_contrast,contrasts(c).conplus);
        % average any multiples
        meanpos = nanmean(groupres.(targetfield)(posind,:,:),1);
        if isempty(contrasts(c).conminus)
            % a vs 0 contrast
            meanneg = zeros(size(meanpos));
        else
            % a vs b contrast
            negind = findconind(groupres.rows_contrast,contrasts(c).conminus);
            assert(~any(intersect(posind,negind)),...
                'same predictor on both sides of contrast');
            meanneg = nanmean(groupres.(targetfield)(negind,:,:),1);
        end
        groupres.(targetfield)(end+1,:,:) = meanpos-meanneg;
        groupres.rows_contrast{end+1} = contrasts(c).name;
        tail{end+1} = contrasts(c).tail;
    end
end

% remove undersized ROIs
n = sum(~isnan(groupres.(targetfield)),3);
% conditions that are all NaN across ROIs
badcon = all(n==0,2);
assert(~all(badcon),'no valid data in targetfield:%s',targetfield);

% for any conditions that are not all NaN, we need the sample size to be
% matched across conditions
firstvalidcon = find(~badcon,1,'first');
r = arrayfun(@(x)all(n(x,:)==n(firstvalidcon,:)),[1:size(n,1)]');
assert(all(r(~badcon)),...
    'mismatched sample size across contrasts for a given ROI');
badroi = n(firstvalidcon,:) < minn;
if any(badroi)
    logstr('removing %d/%d rois with sample size <%d\n',sum(badroi),...
        nroi,minn);
end
if any(badcon)
    logstr('removing %d/%d conditions with no data\n',sum(badcon),ncon);
end

% drop those rois and cons
groupres.cols_roi(badroi) = [];
groupres.rows_contrast(badcon) = [];
tail(badcon) = [];
for t = targets(:)'
    tstr = t{1};
    % this is needed to deal with tail (n by 1 array)
    if size(groupres.(tstr),2) > 1
        groupres.(tstr)(:,badroi,:) = [];
    end
    % this is needed to deal with nfeatures (1 by n array)
    if size(groupres.(tstr),1) > 1
        groupres.(tstr)(badcon,:,:) = [];
    end
end
if isfield(groupres,'custom')
    for c = 1:ncustom
        groupres.custom{c} = indexdim(groupres.custom{c},~badroi,2);
    end
end
nroi = size(groupres.(targetfield),2);
ncon = numel(groupres.rows_contrast);
groupres.z_subject = names;

% now we can get the mean struct 
meanres = struct('rows_contrast',{groupres.rows_contrast},...
    'cols_roi',{groupres.cols_roi},'z_subject',{names});
meanres.tail = tail;
meanres.mean = nanmean(groupres.(targetfield),3);
meanres.n = sum(~isnan(groupres.(targetfield)),3);
meanres.std = nanstd(groupres.(targetfield),[],3);
meanres.sterr = meanres.std ./ sqrt(meanres.n);
% nb r is now the group r stat (on the r stats), and ppara / pperm
% are computed on the subjects too
meanres.t = meanres.mean ./ meanres.sterr;

% p values with correct tail
meanres.ppara = NaN(size(meanres.t));
tind = strcmp(meanres.tail,'right');
meanres.ppara(tind,:) = tcdf(-meanres.t(tind,:),meanres.n(tind,:)-1);
tind = strcmp(meanres.tail,'left');
meanres.ppara(tind,:) = tcdf(meanres.t(tind,:),meanres.n(tind,:)-1);
tind = strcmp(meanres.tail,'both');
% 2 tailed
meanres.ppara(tind,:) = tcdf(-abs(meanres.t(tind,:)),...
    meanres.n(tind,:)-1) * 2;

if isfield(groupres,'nfeatures')
    meanres.nfeatures = nanmean(groupres.nfeatures,3);
    meanres.nfeatures_std = nanstd(groupres.nfeatures,[],3);
end

switch multpdim
    case 'roi'
        % iterate over cons, correct over all ROIs
        fwedim = 1;
        fwelim = ncon;
        tail = meanres.tail;
        threshshape = [ncon 1];
    case 'con'
        % iverate over ROIs, correct over all cons
        fwedim = 2;
        fwelim = nroi;
        % this mode is useful for e.g. RDM multiple comparisons correction,
        % but it does come with a few caveat.
        assert(isequal(meanres.tail{1},meanres.tail{:}),...
            'all cons must have the same tail if multpdim=roi');
        assert(~assumeregister,['multpdim=roi and assumeregister ' ...
            'cannot be combined']);
        threshshape = [1 nroi];
    otherwise
        error('unknown multpdim: %s',multpdim);
end

varsmooth = false;
if ~isempty(varsmoothmask)
    varsmooth = true;
    assert(assumeregister,'must assumeregister if variance smoothing');
    V = spm_vol(varsmoothmask);
    maskxyz = spm_read_vols(V)>0;
    voxsize = vox2mm(V);
    meanres.pseudot = NaN(size(meanres.t));
    meanres.pfwepermvarsm = NaN(size(meanres.t));
end


meanres.pperm = NaN(size(meanres.t));
meanres.pfweperm = NaN(size(meanres.t));
meanres.pfdrthresh = NaN(threshshape);
meanres.medianboot = NaN(size(meanres.t));
meanres.medianste = NaN(size(meanres.t));

% segregate the nulldist/bootdist because these tend to blow up in size
% Key problem: these are often too big to work with. So now we only
% generate them if you asked for them in output
bigvars = false;
if nargout >1 || ~assumeregister
    nulldist = struct('cols_roi',{meanres.cols_roi},'rows_contrast',...
        {meanres.rows_contrast},targetfield,NaN([ncon nroi nperm]));
    bootdist = struct('cols_roi',{meanres.cols_roi},'rows_contrast',...
        {meanres.rows_contrast},targetfield,NaN([ncon nroi nboot]));
    bigvars = true;
end

% for ease of indexing, data is now
% subject by contrast by roi
roidata = shiftdim(groupres.(targetfield),2);
if nperm > 1 || nboot > 0
    % we only have to bother with all this fuss if we are doing
    % some kind of resampling
    if assumeregister
        % life is easy, one model and test per contrast. This is generally
        % faster than looping over rois because we are often dealing
        % with a searchlight volume here (many times more features than
        % contrasts)
        for c = 1:ncon
            fprintf('running %s (%d of %d)...',meanres.rows_contrast{c},...
                c,ncon);
            tstart = clock;

            design = ones([meanres.n(c,1),1],class(roidata));
            assert(numel(unique(meanres.n(c,:)))==1,...
                'n must be the same for all ROIs to use assumeregister');
            % avoid squeeze to prevent weird behaviours
            model = GLM(design,reshape(roidata(:,c,:),[size(roidata,1) ...
                size(roidata,3)]));

            % permutation test
            if nperm > 1
                if varsmooth
                    % variance smooth permutation test
                    [tempnull,tempvar] = permflipsamples(model,nperm,...
                        {'contrast','standarderror'},{[1 nroi],[1 nroi]},1);
                    [meanres.pfwepermvarsm(c,:),meanres.pseudot(c,:)] = ...
                        permpfwe_varsmooth(squeeze(tempnull),...
                        squeeze(tempvar),...
                        maskxyz,'fwhm',varsmoothfwhm,'voxsize',...
                        voxsize,'tail',meanres.tail{c});
                else
                    % use T to avoid max stat being dominated by
                    % high-variance searchlights
                    tempnull = permflipsamples(model,nperm,'tmap',...
                        [1 nroi],1);
                end
                meanres.pperm(c,:) = permpvalue(tempnull,...
                    meanres.tail{c});
                meanres.pfweperm(c,:) = permpfwe(squeeze(tempnull),...
                    meanres.tail{c});
                if bigvars
                    nulldist.(targetfield)(c,:,:) = tempnull;
                end
            end
            try
                meanres.pfdrthresh(c) = FDRthreshold(meanres.ppara(c,:));
            catch
                err = lasterror;
                if ~strcmp(err.identifier,'MATLAB:badsubscript')
                    rethrow(err);
                end
                warning('FDR estimation failed - not adequate p distro?');
            end

            % bootstrap
            if nboot > 0
                tempboot = bootstrapsamples(model,nboot,'fit',[1 nroi]);
                [meanres.medianboot(c,:),meanres.medianste(c,:)] = ...
                    bootprctile(tempboot);
                if model.nsamples < 4
                    fprintf('n<4, rescoring medianste to NaN'\n');
                    meanres.medianste(c,:) = NaN;
                end
                if bigvars
                    bootdist.(targetfield)(c,:,:)  = tempboot;
                end
            end

            fprintf(' finished in %s.\n',seconds2str(etime(clock,tstart)));
        end
    else
        % have to do each ROI separately to get correct n. Instead we put
        % all the contrasts in feature dimension for the model for some
        % speed up compared to fully nested loop case.
        for r = 1:nroi
            % indexing the first row works because n is the same for
            % each contrast
            design = ones(meanres.n(1,r),1);
            thisroidata = roidata(:,:,r);
            % any and all are assumed to be equivalent due to test
            % above (nb, magically design is already the right length)
            goodroi = ~any(isnan(thisroidata),2);
            model = GLM(design,thisroidata(goodroi,:));

            % permutation test
            if nperm > 1
                % nb permute T, not raw estimates, to avoid issues with
                % variance inhomogeneities across roi/con.
                nulldist.(targetfield)(:,r,:) = permflipsamples(model,...
                    nperm,'tmap',[1 ncon],1);
            end

            % bootstrap
            if nboot > 0
                bootdist.(targetfield)(:,r,:) = bootstrapsamples(model,...
                    nboot,'fit',[1 ncon]);
            end
        end
        for c = 1:fwelim
            % get the current ROI or condition depending on multpdim
            outind = dimindices(size(meanres.t),c,fwedim);
            if nperm > 1
                innull = indexdim(nulldist.(targetfield),c,fwedim);
                meanres.pperm(outind) = permpvalue(innull,tail{c});
                % unlike pperm, pfweperm needs matrix input
                % so remove the fwedim, while avoiding squeeze
                sz = size(innull);
                sz(fwedim) = [];
                innull = reshape(innull,sz);
                meanres.pfweperm(outind) = permpfwe(innull,...
                    tail{c});
            end
            % don't need tail here because the ppara field is already
            % tailed
            try
                meanres.pfdrthresh(c) = FDRthreshold(indexdim(...
                    meanres.ppara,c,fwedim));
            catch
                err = lasterror;
                if ~strcmp(err.identifier,'MATLAB:badsubscript')
                    rethrow(err);
                end
                warning('FDR estimation failed - not adequate p distro?');
            end
        end
        if nboot > 0
            [meanres.medianboot,meanres.medianste] = bootprctile(...
                bootdist.(targetfield));
            if model.nsamples < 4
                fprintf('n<4, rescoring medianste to NaN'\n');
                meanres.medianste(:) = NaN;
            end
        end
    end
end

% do the customs as well
for nc = 1:ncustom_rfx
    for r = 1:nroi
        meanres.custom{nc,r} = feval(customfun{nc},roidata(:,:,r),...
            meanres.rows_contrast);
    end
end

if keepnans && any(badcon)
    logstr('returning nan conditions to output\n');
    meanres = upcastdata(meanres,~badcon,targetfield);
    groupres = upcastdata(groupres,~badcon,targetfield);
end

%% SUBFUNCTIONS
function meanres = upcastdata(meanres,validind,targetfield)

outsize = size(meanres.(targetfield));
nvalid = size(validind,1);
for fn = fieldnames(meanres)'
    fnstr = fn{1};
    thissize = size(meanres.(fnstr));
    if isequal(outsize,thissize)
        % need to upcast
        newv = NaN([nvalid outsize(2:end)]);
        newv(repmat(validind,[1 outsize(2:end)])) = ...
            meanres.(fnstr);
        meanres.(fnstr) = newv;
    elseif outsize(1)==thissize(1) && thissize(2)==1
        % vector
        newv = cell(nvalid,1);
        newv(validind) = meanres.(fnstr);
        meanres.(fnstr) = newv;
    else
        logstr('cannot parse %s, skipping...\n',fnstr);
    end
end
