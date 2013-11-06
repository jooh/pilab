% Perform information-based T mapping on each of the ROIs in the rois
% volume instance, using discriminants in contrasts fitted using designvol
% and epivol. 
%
% Named inputs:
%
% glmclass: char defining CovGLM sub-class (e.g. 'CovGLM')
% glmvarargs: any additional arguments for GLM (e.g. k for RidgeGLM)
%
% (used here)
% split: [] passed to cvgroups to define crossvalidation
% nperm: 1 number of permutations to run (1 being the first, true order)
% nboot: 0 number of bootstrap resamples of the true model
%
% [res,nulldist,bootdist] = roidata_lindisc(rois,designvol,epivol,contrasts,varargin)
function [res,nulldist,bootdist] = roidata_lindisc(rois,designvol,epivol,contrasts,varargin)

ts = varargs2structfields(varargin,struct(...
    'split',[],'glmclass','GLM','glmvarargs',[],...
    'nperm',1,'nboot',0));

if ~iscell(ts.split)
    % so we can easily deal to cvgroup field later
    ts.split = num2cell(ts.split);
end

% we now assume that the inputs are super clean (so that we can assume that
% the data result will be in register with rois)
assert(~any(isnan(epivol.data(:))),'nans in epivol');
assert(isequal(epivol.mask,rois.mask),...
    'mismatched masks in epivol and rois');

% set resampling once so all null distributions are yoked and can be
% used to obtain difference distributions between ROIs or classifiers
perminds = permuteindices(designvol.desc.samples.nunique.chunks,ts.nperm);
bootinds = bootindices(designvol.desc.samples.nunique.chunks,ts.nboot);

ncon = numel(contrasts);
basedat = NaN([ncon rois.nsamples]);
cols_roi = {asrow(rois.meta.samples.names)};
rows_contrast = {ascol({contrasts.name})};

res = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    't',basedat,'pperm',basedat,'medianboot',basedat,...
    'medianste',basedat,'ppara',basedat,'nfeatures',...
    NaN([1,rois.nsamples]));

% segregate the nulldist/bootdist because these tend to blow up in size
nulldist = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    't',NaN([ncon rois.nsamples ts.nperm]));
bootdist = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    't',NaN([ncon rois.nsamples ts.nboot]));

% can't parfor here because epivol is likely too big to pass around
for r = 1:rois.nsamples
    % skip empty rois (these come out as NaN)
    validvox = full(rois.data(r,:)~=0);
    if ~any(validvox)
        % empty roi
        fprintf('no valid voxels, skipping..\n');
        continue
    end

    % make the GLM instance for this ROI
    % (we will need to strip out 0 cons eventually for speed - maybe move
    % this inside contrast loop and apply a skipcon that gets applied just
    % before volume stage?)
    model = vol2glm(designvol,epivol(:,validvox),ts.glmclass,...
        ts.glmvarargs{:});
    % nb model can still have multiple array entries
    [model.cvgroup] = ts.split{:};
    res.nfeatures(r) = model.nfeatures;

    % permute design matrix
    if ts.nperm > 1
        fprintf('running %d permutations\n',ts.nperm);
    end

    % get discriminants
    for c = 1:ncon
        % this extra if/else is a little redundant but it's good to avoid
        % parfor overhead whenever possible
        if ts.nperm < 2
            raw = cvcrossclassificationrun(model,'discriminant',...
                'infotmap',[],{contrasts(c).traincon},...
                {contrasts(c).testcon});
            % average splits and discriminants
            nulldist.t(c,r,1) = mean(raw(:));
        else
            % permutation test
            ptemp = nulldist.t(c,r,:);
            parfor p = 1:ts.nperm
                % draw a model (or just the original if p==1)
                permmodel = drawpermrun(model,perminds(p,:));
                raw = cvcrossclassificationrun(permmodel,'discriminant',...
                    'infotmap',[],{contrasts(c).traincon},...
                    {contrasts(c).testcon});
                % average splits and discriminants
                ptemp(1,1,p) = mean(raw(:));
            end % p ts.nperm
            nulldist.t(c,r,:) = ptemp;
        end

        if ts.nboot > 0
            % bootstrap
            btemp = bootdist.t(c,r,:);
            parfor b = 1:ts.nboot
                bootmodel = model(bootinds(b,:));
                raw = cvcrossclassificationrun(bootmodel,'discriminant',...
                    'infotmap',[],{contrasts(c).traincon},...
                    {contrasts(c).testcon});
                % average splits and discriminants
                btemp(1,1,b) = mean(raw(:));
            end
            bootdist.t(c,r,:) = btemp;
        end
    end % c ncon
end % r rois.nsamples

% get result for true model
res.t = nulldist.t(:,:,1);
if ts.nperm > 1
    res.pperm = permpvalue(nulldist.t);
end
% one tailed p values (nb this calculation isn't exactly right since we are
% actually working with an 'average' t across many splits and
% discriminants. However, I think the rationale is that p is actually
% conservative compared to what it should be. I guess comparisons with the
% perm p values will show.
res.ppara = tcdf(-res.t,sum(vertcat(model.nsamples))-1);
[res.medianboot,res.medianste] = bootprctile(bootdist.t);
