% Perform rsa mapping on each of the RDMs in the disvol compared to each
% model RDMs in the predictors struct.
%
% INPUTS
% disvol        BaseVolume instance containing vectorised dissimilarities
%                   (or cell array for multiple run support)
% predictors    model RDMs in struct array form
%
% NAMED INPUTS
%
% ARGUMENT          DEFAULT DESCRIPTION
% rsaclass          RankRSA RSA subclass to use (default RankRSA)
% rsaclassargs      {}      any varargins for rsaclass
% nperm             1       number of permutations to run (1 being the
%                               first, true order)
% nboot             0       number of bootstrap resamples of the true model
% splitrsaclass     ''      RSA subclass to use for split data RDMs
% splitrsaclassargs {}      any varargins for splitrsaclass
% fitmethod         'fit'   RSA class method to call to obtain result
% fitmethodargs     {}      any additional arguments for fitmethod
% customfun         []      catch-all custom processing. Gets called as
%                               res.custom{c} = feval(customfun,model);
%
% [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,[varargin])
function [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,varargin)

getArgs(varargin,{'rsaclass','RankRSA','nperm',1,'nboot',0,...
    'rsaclassargs',{},'splitrsaclass','',...
    'splitrsaclassargs',{},'fitmethod','fit','fitmethodargs',{},...
    'customfun',[]});

if ~iscell(rsaclassargs)
    % this extra bit of flexibility helps with cases where you just want to
    % set a single extra parameter without needing to wrap it in curlies
    rsaclassargs = {rsaclassargs};
end

if ~iscell(splitrsaclassargs)
    splitrsaclassargs = {splitrsaclassargs};
end

if ~iscell(disvol)
    disvol = {disvol};
end

% we now assume that the inputs are super clean (so that we can assume that
% the data result will be in register with rois)
assert(~any(isnan(disvol{1}.data(:))),'nans in disvol');

npredictor = numel(predictors);
basedat = NaN([npredictor disvol{1}.nfeatures]);
cols_roi = {asrow(disvol{1}.meta.features.names)};
rows_contrast = {ascol({predictors.name})};

res = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',basedat,'pperm',basedat,'medianboot',basedat,...
    'medianste',basedat);

if isfield(disvol{1}.meta.features,'nfeatures')
    res.nfeatures = disvol{1}.meta.features.nfeatures;
end

% segregate the nulldist/bootdist because these tend to blow up in size
nulldist = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',NaN([npredictor disvol{1}.nfeatures nperm]));
bootdist = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',NaN([npredictor disvol{1}.nfeatures nboot]));

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
    model = cellfun(@(thisdis)feval(constructor,predictors(c).RDM,...
        thisdis.data,args{:}),disvol,'uniformoutput',0);
    if iscell(model)
        model = cat(1,model{:});
    end

    res.r(c,:) = feval(fitmethod,model,fitmethodargs{:});

    if nperm > 1
        nulldist.r(c,:,:) = permutesamples(model,nperm,fitmethod,[],...
            fitmethodargs{:});
    end

    if nboot > 0
        bootdist.r(c,:,:) = bootstrapsamples(model,nboot,fitmethod,[],...
            fitmethodargs{:});
    end

    if ~isempty(customfun)
        res.custom{c} = feval(customfun,model);
    end

end % c npredictor

if nperm > 1
    res.pperm = permpvalue(nulldist.t);
end
if nboot > 0
    [res.medianboot,res.medianste] = bootprctile(bootdist.t);
end
