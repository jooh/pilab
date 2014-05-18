% do a random effects analysis on the data in struct array subres (one
% entry per subject). Makes various assumptions about how subres is
% organised - see roidata_rsa and roidata_lindisc for reference.
%
% Named inputs:
% nperm: 1 
% nboot: 0
% targetfield: 't'
% transfun: string that gets feval'ed or a function handle to apply some
% transform to the data before doing the rfx analysis (e.g. 'atanh' for
% fisher transform, or something like @(x)x-.5 to mean correct proportion
% correct classification performance).
% assumeregister: false (if true, we just stack the single subject results,
% consequences be damned)
% varsmoothmask: []
% varsmoothfwhm: 8 (mm)
% contrasts: calculate differences between particular conditions. char that
% gets evaled or struct with fields:
%   name
%   conplus
%   conminus
%   tail
%
% [meanres,nulldist,bootdist] = roidata_rfx(subres,[varargin])
function [meanres,nulldist,bootdist] = roidata_rfx(subres,varargin)
    
getArgs(varargin,{'nperm',1,'nboot',0,'targetfield','t','transfun',[],...
    'assumeregister',false,'varsmoothmask',[],'varsmoothfwhm',8,...
    'contrasts',struct('name',{},'conplus',{},'conminus',{},'tail',{})});

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

if ~isempty(contrasts)
    if ischar(contrasts)
        contrasts = feval(contrasts);
    end
    assert(isstruct(contrasts),'contrasts must be a struct');
end

% but we do assume that contrasts are identical.  we could support
% this case in theory but let's not for now since it likely
% indicates a user error.
concell = {subres.rows_contrast};
assert(isequal(concell{:}),...
    'different predictors in different subjects');
concell = concell{1};
ncon = length(concell);
nsub = length(subres);
dat = NaN([ncon nroi nsub]);
groupres = struct('rows_contrast',{subres(1).rows_contrast},...
    'cols_roi',{urois},targetfield,dat,'nfeatures',NaN([1 nroi nsub]));

% populate the groupres struct
if assumeregister
    % hey, this is easy
    groupres.(targetfield) = cat(3,subres.(targetfield));
    assert(isreal(groupres.(targetfield)),'imaginary output!');
    assert(~any(isnan(groupres.(targetfield)(:))),...
        'nans in targetfield not supported in assumeregister mode');
else
    % have to be careful not to average over different rois
    for s = 1:nsub
        for r = 1:length(subres(s).cols_roi)
            thisroi = subres(s).cols_roi{r};
            indroi = strcmp(groupres.cols_roi,thisroi);
            groupres.(targetfield)(:,indroi,s) = ...
                subres(s).(targetfield)(:,r);
            if isfield(groupres,'nfeatures')
                groupres.nfeatures(1,indroi,s) = subres(s).nfeatures(1,r);
            end
        end
    end
end

ptails = repmat({'right'},[ncon,1]);

% add in the contrasts here
if ~isempty(contrasts)
    for c = 1:length(contrasts)
        % find the rows corresponding to each condition
        [~,posind] = intersect(groupres.rows_contrast,...
            contrasts(c).conplus);
        if isempty(posind)
            [~,posind] = intersect(groupres.rows_contrast,cellfun(...
                @(thiscon)['rsa_r_' thiscon],contrasts(c).conplus,...
                'uniformoutput',false));
        end
        assert(~isempty(posind),'did not find %s',contrasts(c).conplus);
        [~,negind] = intersect(groupres.rows_contrast,...
            contrasts(c).conminus);
        if isempty(negind)
            [~,negind] = intersect(groupres.rows_contrast,cellfun(...
                @(thiscon)['rsa_r_' thiscon],contrasts(c).conminus,...
                'uniformoutput',false));
        end
        assert(~isempty(negind),'did not find %s',contrasts(c).conminus);
        assert(~any(intersect(posind,negind)),...
            'same predictor on both sides of contrast');
        % average any multiples
        meanpos = nanmean(groupres.(targetfield)(posind,:,:),1);
        meanneg = nanmean(groupres.(targetfield)(negind,:,:),1);
        groupres.(targetfield)(end+1,:,:) = meanpos-meanneg;
        groupres.rows_contrast{end+1} = contrasts(c).name;
        ptails{end+1} = contrasts(c).tail;
    end
end

ncon = numel(groupres.rows_contrast);

if ~isempty(transfun)
    if ischar(transfun)
        % apply some transform, e.g. fisher
        groupres.(targetfield) = feval(transfun,groupres.(targetfield));
    else
        % direct call on function handle
        groupres.(targetfield) = transfun(groupres.(targetfield));
    end
end

% now we can get the mean struct 
meanres = struct('rows_contrast',{groupres.rows_contrast},...
    'cols_roi',{groupres.cols_roi},'z_subject',{{subres.name}});
meanres.ptails = ptails;
meanres.mean = nanmean(groupres.(targetfield),3);
meanres.n = sum(~isnan(groupres.(targetfield)),3);
meanres.std = nanstd(groupres.(targetfield),[],3);
meanres.sterr = meanres.std ./ sqrt(meanres.n);
% nb r is now the group r stat (on the r stats), and ppara / pperm
% are computed on the subjects too
meanres.t = meanres.mean ./ meanres.sterr;

% p values with correct tail
meanres.ppara = NaN(size(meanres.t));
tind = strcmp(meanres.ptails,'right');
meanres.ppara(tind,:) = tcdf(-meanres.t(tind,:),meanres.n(tind,:)-1);
tind = strcmp(meanres.ptails,'left');
meanres.ppara(tind,:) = tcdf(meanres.t(tind,:),meanres.n(tind,:)-1);
tind = strcmp(meanres.ptails,'both');
% 2 tailed
meanres.ppara(tind,:) = tcdf(-abs(meanres.t(tind,:)),...
    meanres.n(tind,:)-1) * 2;

if isfield(groupres,'nfeatures')
    meanres.nfeatures = nanmean(groupres.nfeatures,3);
    meanres.nfeatures_std = nanstd(groupres.nfeatures,[],3);
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

% check that sample size is the same for all contrasts by comparing
% each to the first
r = bsxfun(@eq,meanres.n,meanres.n(1,:));
assert(all(r(:)),...
    'mismatched sample size across contrasts for a given ROI');

meanres.pperm = NaN(size(meanres.t));
meanres.pfweperm = NaN(size(meanres.t));
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
            model = GLM(design,squeeze(roidata(:,c,:)));

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
                        voxsize,'tail',meanres.ptails{c});
                else
                    % use T to avoid max stat being dominated by
                    % high-variance searchlights
                    tempnull = permflipsamples(model,nperm,'tmap',...
                        [1 nroi],1);
                end
                meanres.pperm(c,:) = permpvalue(tempnull,...
                    meanres.ptails{c});
                meanres.pfweperm(c,:) = permpfwe(squeeze(tempnull)',...
                    meanres.ptails{c});
                if bigvars
                    nulldist.(targetfield)(c,:,:) = tempnull;
                end
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
                nulldist.(targetfield)(:,r,:) = permflipsamples(model,...
                    nperm,'fit',[1 ncon]);
            end

            % bootstrap
            if nboot > 0
                bootdist.(targetfield)(:,r,:) = bootstrapsamples(model,...
                    nboot,'fit',[1 ncon]);
            end
        end
        % nb, when ~assumeregister we basically assume bigvars will be ok
        % (typically this is an ROI analysis with manageable data size)
        if nperm > 1
            for c = 1:ncon
                meanres.pperm(c,:) = permpvalue(...
                    nulldist.(targetfield)(c,:,:),meanres.ptails{c});
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

