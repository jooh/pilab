% Perform rsa mapping on each of the RDMs in the disvol compared to each
% model RDMs in the predictors struct.
%
% INPUTS
% disvol        Volume instance containing vectorised dissimilarities
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
% bootmeth          bootstrapsamples    method to call for bootstrapping
% bootprep          preparesampleboots  method to initialise bootinds
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
%
% [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,[varargin])
function [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,varargin)

getArgs(varargin,{'rsaclass','RankRSA','nperm',1,'nboot',0,...
    'rsaclassargs',{},'splitrsaclass','',...
    'splitrsaclassargs',{},'fitmethod','fit','fitmethodargs',{},...
    'bootmeth','bootstrapsamples','bootprep','preparesampleboots',...
    'customfun',[],...
    'contrasts',emptystruct('name','conplus','conminus','tail'),...
    'customfits',emptystruct('name','funhand','cons','tail')});

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

% we now assume that the inputs have no missing ROIs (to preserve the
% registration across subjects)
assert(~any(all(isnan(disvol{1}.data),1)),'nan ROIs in disvol');

npredictor = numel(predictors);
basedat = NaN([npredictor disvol{1}.nfeatures],class(disvol{1}.data));
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
    'r',NaN([npredictor disvol{1}.nfeatures nperm],class(disvol{1}.data)));
bootdist = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',NaN([npredictor disvol{1}.nfeatures nboot],class(disvol{1}.data)));

permind = [];
bootind = [];
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
        if isempty(permind)
            % cache the sample perms so that all conditions are
            % permuted similarly (enabled contrasts between null
            % distributions later).
            permind = preparesampleperms(model,nperm);
        end
        nulldist.r(c,:,:) = permutesamples(model,permind,fitmethod,[],...
            fitmethodargs{:});
    end

    if nboot > 0
        if isempty(bootind)
            bootind = feval(bootprep,model,nboot);
        end
        bootdist.r(c,:,:) = feval(bootmeth,model,bootind,fitmethod,[],...
            fitmethodargs{:});
    end

    if ~isempty(customfun)
        res.custom{c} = feval(customfun,model);
    end

end % c npredictor

if ~isempty(customfits)
    error('customfits input is not yet supported!')
end

% add in the contrasts
if ~isempty(contrasts)
    if ischar(contrasts)
        contrasts = feval(contrasts,res.rows_contrast);
    end
    assert(isstruct(contrasts),'contrasts must be a struct');
    for c = 1:length(contrasts)
        % find the rows corresponding to each condition
        posind = findconind(res.rows_contrast,contrasts(c).conplus);
        % average any multiples
        meanpos = nanmean(res.r(posind,:),1);
        if isempty(contrasts(c).conminus)
            % a vs 0 contrast
            meanneg = zeros(size(meanpos));
        else
            % a vs b contrast
            negind = findconind(res.rows_contrast,contrasts(c).conminus);
            assert(~any(intersect(posind,negind)),...
                'same predictor on both sides of contrast');
            meanneg = nanmean(res.r(negind,:),1);
        end
        res.r(end+1,:) = meanpos-meanneg;
        res.rows_contrast{end+1} = contrasts(c).name;
        if nperm > 1
            % permutation test for contrast by averaging permutation stats
            if isempty(contrasts(c).conminus)
                % just average
                nulldist.r(end+1,:,:) = nanmean(...
                    nulldist.r(posind,:,:),1);
            else
                nulldist.r(end+1,:,:) = nanmean(...
                    nulldist.r(posind,:,:),1) - ...
                    nanmean(nulldist.r(negind,:,:),1);
            end
            % this is a bit hacky but it's more convenient to flip the null
            % distribution around than to feed lots of custom arguments to
            % permpvalue below.
            switch contrasts(c).tail
                case 'right'
                    % no problem
                case 'left'
                    nulldist.r(end,:,:) = nulldist.r(end,:,:) * -1;
                case 'both'
                    nulldist.r(end,:,:) = abs(nulldist.r(end,:,:));
                otherwise
                    error('unknown tail argument: %s',contrasts(c).tail);
            end
        end % if nperm > 1
        if nboot > 1
            % bootstrap distribution for contrast
            if isempty(contrasts(c).conminus)
                % just average
                bootdist.r(end+1,:,:) = nanmean(...
                    bootdist.r(posind,:,:),1);
            else
                bootdist.r(end+1,:,:) = nanmean(...
                    bootdist.r(posind,:,:),1) - ...
                    nanmean(bootdist.r(negind,:,:),1);
            end
        end
    end
end

if nperm > 1
    res.pperm = permpvalue(nulldist.r);
end
if nboot > 0
    [res.medianboot,res.medianste] = bootprctile(bootdist.r);
end
