% Perform rsa mapping on each of the RDMs in the disvol compared to each
% model RDMs in the predictors struct.
%
% Named inputs (default and description):
%
% rsaclass: RankRSA RSA subclass to use (default RankRSA)
% rsaclassargs: {} any varargins for rsaclass
% nperm: 1 number of permutations to run (1 being the first, true order)
% nboot: 0 number of bootstrap resamples of the true model
% splitrsaclass: Custom RSA subclass to use for split data RDMs (default
% SplitRankRSA)
% splitrsaclassargs: {} any varargins for splitrsaclass
%
% [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,[varargin])
function [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,varargin)

getArgs(varargin,{'rsaclass','RankRSA','nperm',1,'nboot',0,...
    'rsaclassargs',{},'splitrsaclass','SplitRankRSA','splitrsaclassargs',{}});

if ~iscell(rsaclassargs)
    % this extra bit of flexibility helps with cases where you just want to
    % set a single extra parameter without needing to wrap it in curlies
    rsaclassargs = {rsaclassargs};
end

if ~iscell(splitrsaclassargs)
    splitrsaclassargs = {splitrsaclassargs};
end

% we now assume that the inputs are super clean (so that we can assume that
% the data result will be in register with rois)
assert(~any(isnan(disvol.data(:))),'nans in disvol');

npredictor = numel(predictors);
basedat = NaN([npredictor disvol.nfeatures]);
cols_roi = {asrow(disvol.meta.features.names)};
rows_contrast = {ascol({predictors.name})};

res = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',basedat,'pperm',basedat,'medianboot',basedat,...
    'medianste',basedat);

if isfield(disvol.meta.features,'nfeatures')
    res.nfeatures = disvol.meta.features.nfeatures;
end

% segregate the nulldist/bootdist because these tend to blow up in size
nulldist = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',NaN([npredictor disvol.nfeatures nperm]));
bootdist = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',NaN([npredictor disvol.nfeatures nboot]));

% is this stored somewhere in disvol? Suspect it is.
% res.nfeatures(r) = model.nfeatures;
for c = 1:npredictor
    if issplitdatardm(predictors(c).RDM)
        % split class instance
        constructor = splitrsaclass;
        args = splitrsaclassargs;
    else
        % regular class instance
        constructor = rsaclass;
        args = rsaclassargs;
    end
    % setup the instance
    model = feval(constructor,predictors(c).RDM,disvol.data,args{:});

    res.r(c,:) = fit(model);

    if nperm > 1
        nulldist.r(c,:,:) = permutesamples(model,nperm,'fit');
    end

    if nboot > 0
        bootdist.r(c,:,:) = bootstrapsamples(model,nboot,'fit');
    end

end % c npredictor

if nperm > 1
    res.pperm = permpvalue(nulldist.t);
end
if nboot > 0
    [res.medianboot,res.medianste] = bootprctile(bootdist.t);
end
