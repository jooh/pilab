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
% permdis           []      ndis by nroi by nperm matrix of permuted RDMs
%                               from e.g. GLM.permuteruns. We use these to
%                               build nulldist in place of the RSA
%                               class-based permutation methods. If this
%                               argument is available nperm must be < 2.
% nperm             1       number of permutations to run (1 being the
%                               first, true order)
% nboot             0       number of bootstrap resamples of the true model
% splitrsaclass     ''      RSA subclass to use for split data RDMs
% splitrsaclassargs {}      any varargins for splitrsaclass
% fitmethod         'fit'   RSA class method to call to obtain result
% fitmethodargs     {}      any additional arguments for fitmethod
% customfun         []      catch-all custom processing. Gets called as
%                               res.custom{c,n} = feval(customfun{n},model);
% bootmeth          bootstrapsamples    method to call for bootstrapping
% bootprep          preparesampleboots  method to initialise bootinds
% defaulttail       'right' tail for tests
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
%       res.r(end+1,:,:) = customfits(c).funhand(...
%       res.r(conind,:,:),customfits(c));
%
% [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,[varargin])
function [res,nulldist,bootdist] = roidata_rsa(disvol,predictors,varargin)

getArgs(varargin,{'rsaclass','RankRSA','nperm',1,'nboot',0,...
    'rsaclassargs',{},'splitrsaclass','',...
    'splitrsaclassargs',{},'fitmethod','fit','fitmethodargs',{},...
    'bootmeth','bootstrapsamples','bootprep','preparesampleboots',...
    'customfun',[],'defaulttail','right',...
    'contrasts',emptystruct('name','conplus','conminus','tail'),...
    'customfits',emptystruct('name','funhand','cons','tail'),...
    'permdis',[]});

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

if ~iscell(customfun)
    if isempty(customfun)
        customfun = {};
    else
        customfun = {customfun};
    end
end
ncustom = numel(customfun);

npermdis = 1;
if ~isempty(permdis)
    if ~iscell(permdis)
        permdis = {permdis};
    end
    permdissz = size(permdis{1});
    if ndims(permdis{1}) > 2
        assert(nperm < 2,'nperm must be <2 if permdis is provided');
        npermdis = permdissz(3);
        assert(isequal(permdissz(1:2),...
            [disvol{1}.nsamples disvol{1}.nfeatures]),['size mismatch: ' ...
            'disvol{1} - permdis{1}']);
        assert(numel(permdis) == numel(disvol),['numel mismatch: disvol - '...
            'permdis']);
    end
end

% we now assume that the inputs have no missing ROIs (to preserve the
% registration across subjects)
assert(~any(all(isnan(disvol{1}.data),1)),'nan ROIs in disvol');

npredictor = numel(predictors);
basedat = NaN([npredictor disvol{1}.nfeatures],class(disvol{1}.data));

cols_roi = {asrow(disvol{1}.meta.features.names)};
if numel(cols_roi)==1 && isempty(cols_roi{1})
    cols_roi = {cell(1,disvol{1}.nfeatures)};
end
rows_contrast = {ascol({predictors.name})};

dismat = cellfun(@(thisdis)thisdis.data,disvol,'uniformoutput',0);

res = struct('cols_roi',cols_roi,'rows_contrast',rows_contrast,...
    'r',basedat);

if isfield(disvol{1}.meta.features,'nfeatures')
    res.nfeatures = disvol{1}.meta.features.nfeatures;
end

% segregate the nulldist/bootdist because these tend to blow up in size
nulldist = struct(...
    'r',NaN([npredictor disvol{1}.nfeatures max([nperm npermdis])],...
    class(disvol{1}.data)));
bootdist = struct(...
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

    % get the (true) estimates
    allarg = {predictors(c).RDM,constructor,args,fitmethod,fitmethodargs};
    [res.r(c,:),model] = constructandfit(dismat,1,allarg{:});

    if npermdis > 1
        % RDM-based permutation test
        % temp assignment to keep parfor happy
        r = cell(1,npermdis);
        parfor p = 1:npermdis
            % NB the p gives the z index into each mat contained in the
            % cell array dismat
            %nulldist.r(c,:,p) = constructandfit(permdis,p,allarg{:});
            r{p} = constructandfit(permdis,p,allarg{:});
        end
        nulldist.r(c,:,:) = cat(3,r{:});
    end

    if nperm > 1
        if isempty(permind)
            % cache the sample perms so that all conditions are
            % permuted similarly (enables contrasts between null
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

    for nc = 1:ncustom
        res.custom{c,nc} = feval(customfun{nc},model);
    end

end % c npredictor

res.tail = repmat({defaulttail},[npredictor 1]);

if ~isempty(customfits)
    if ischar(customfits)
        customfits = feval(customfits,res.rows_contrast);
    end
    assert(isstruct(customfits),'customfits must be a struct');
    for c = 1:length(customfits)
        if ischar(customfits(c).funhand)
            customfits(c).funhand = str2func(customfits(c).funhand);
        end
        conind = findconind(res.rows_contrast,customfits(c).cons);
        % fit
        res.r(end+1,:,:) = customfits(c).funhand(...
            res.r(conind,:),customfits(c));
        res.rows_contrast{end+1} = customfits(c).name;
        res.tail{end+1} = customfits(c).tail;
        if max([nperm npermdis]) > 1
            r = cell(1,size(nulldist.r,3));
            parfor n = 1:size(nulldist.r,3)
                r{n} = customfits(c).funhand(...
                    nulldist.r(conind,:,n),customfits(c));
            end
            nulldist.r(end+1,:,:) = cat(3,r{:});
        end % if max([nperm npermdis]) > 1
        if nboot > 1
            r = cell(1,size(bootdist.r,3));
            for b = 1:size(bootdist.r,3)
                r{b} = customfits(c).funhand(...
                    bootdist.r(conind,:,b),customfits(c));
            end
            bootdist.r(end+1,:,:) = cat(3,r{:});
        end
    end
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
        res.tail{end+1} = contrasts(c).tail;
        if max([nperm npermdis]) > 1
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
        end % if max([nperm npermdis]) > 1
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

bootdist.rows_contrast = res.rows_contrast;
bootdist.cols_roi = res.cols_roi;
nulldist.rows_contrast = res.rows_contrast;
nulldist.cols_roi = res.cols_roi;

if max([nperm npermdis]) > 1
    res.pperm = NaN(size(res.r));
    for ttarget = {'right','left','both'}
        tstr = ttarget{1};
        tind = strcmp(res.tail,tstr);
        res.pperm(tind,:,:) = permpvalue(nulldist.r(tind,:,:,:),tstr);
    end
end
if nboot > 0
    [res.medianboot,res.medianste] = bootprctile(bootdist.r);
end

function [r,model] = constructandfit(dismat,zind,predictor,constructor,args,fitmeth,fitmethargs)

% setup the instance
model = cellfun(@(thisdis)feval(constructor,predictor,thisdis(:,:,zind),...
    args{:}),dismat,'uniformoutput',0);
if iscell(model)
    model = cat(1,model{:});
end

r = feval(fitmeth,model,fitmethargs{:});

function res = fitslope(data,customfit)

[nsamp,nroi,nsub] = size(data);
res = NaN([1 nroi nsub]);
for r = 1:nroi
    valid = ~isnan(data(1,r,:));
    roidata = squeeze(data(:,r,valid));
    x = [ones(nsamp,1) customfit.x];
    b = x \ roidata;
    res(1,r,valid) = b(2,:);
end

function res = fitrz(data,customfit)

[nsamp,nroi,nsub] = size(data);
res = NaN([1 nroi nsub]);
for r = 1:nroi
    valid = ~isnan(data(1,r,:));
    roidata = squeeze(data(:,r,valid));
    x = customfit.x(:);
    % fisher-transformed correlation coefficient
    res(1,r,valid) = atanh(zscore(x) \ zscore(roidata));
end
