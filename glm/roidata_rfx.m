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
%
% [meanres,nulldist,bootdist] = roidata_rfx(subres,[varargin])
function [meanres,nulldist,bootdist] = roidata_rfx(subres,varargin)
    
getArgs(varargin,{'nperm',1,'nboot',0,'targetfield','t','transfun',[],...
    'assumeregister',false});

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
    assert(~any(isnan(groupres.(targetfield)(:))),...
        'nans in targetfield not supported in assumeregister mode');
else
    % have to be careful not to average over different rois
    for s = 1:nsub
        for r = 1:length(subres(s).cols_roi)
            thisroi = subres(s).cols_roi{r};
            indroi = strcmp(groupres.cols_roi,thisroi);
            groupres.(targetfield)(:,indroi,s) = subres(s).(targetfield)(:,r);
            if isfield(groupres,'nfeatures')
                groupres.nfeatures(1,indroi,s) = subres(s).nfeatures(1,r);
            end
        end
    end
end

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
meanres.mean = nanmean(groupres.(targetfield),3);
meanres.n = sum(~isnan(groupres.(targetfield)),3);
meanres.std = nanstd(groupres.(targetfield),[],3);
meanres.sterr = meanres.std ./ sqrt(meanres.n);
% nb r is now the group r stat (on the r stats), and ppara / pperm
% are computed on the subjects too
meanres.t = meanres.mean ./ meanres.sterr;
meanres.ppara = tcdf(-meanres.t,meanres.n-1);
if isfield(groupres,'nfeatures')
    meanres.nfeatures = nanmean(groupres.nfeatures,3);
    meanres.nfeatures_std = nanstd(groupres.nfeatures,[],3);
end

% check that sample size is the same for all contrasts by comparing
% each to the first
r = bsxfun(@eq,meanres.n,meanres.n(1,:));
assert(all(r(:)),...
    'mismatched sample size across contrasts for a given ROI');

% now we want to a) bootstrap, b) permutation test. Hm. The
% preferred and maximally future-proof version would be to create a
% GLM with one feature per ROI and run it separately on each
% contrast (with a constant). That way we can use sample-based
% permutation and bootstrap methods. The catch is that we have lots
% of NaNs in the data here so we would probably wind up with a
% different model for each contrast and ROI. This is a pain but it
% is the 'correct' way to do things.
meanres.pperm = NaN(size(meanres.t));
meanres.medianboot = NaN(size(meanres.t));
meanres.medianste = NaN(size(meanres.t));

% segregate the nulldist/bootdist because these tend to blow up in size
nulldist = struct('cols_roi',{meanres.cols_roi},'rows_contrast',...
    {meanres.rows_contrast},targetfield,NaN([ncon nroi nperm]));
bootdist = struct('cols_roi',{meanres.cols_roi},'rows_contrast',...
    {meanres.rows_contrast},targetfield,NaN([ncon nroi nboot]));

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
            design = ones(meanres.n(c,1),1);
            assert(numel(unique(meanres.n(c,:)))==1,...
                'n must be the same for all ROIs to use assumeregister');
            model = GLM(design,squeeze(roidata(:,c,:)));

            % permutation test
            if nperm > 1
                nulldist.(targetfield)(c,:,:) = permflipsamples(model,...
                    nperm,'fit',[1 nroi]);
            end

            % bootstrap
            if nboot > 0
                bootdist.(targetfield)(c,:,:) = bootstrapsamples(model,...
                    nboot,'fit',[1 nroi]);
            end

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
    end
end

if nperm > 1
    meanres.pperm = permpvalue(nulldist.(targetfield));
end
if nboot > 0
    [meanres.medianboot,meanres.medianste] = bootprctile(...
        bootdist.(targetfield));
end
