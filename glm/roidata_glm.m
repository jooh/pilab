% Perform univariate modelling of the input designvol and epivol. Similar
% functionality to roidata_rsa. designvol can be a Volume instance
% containing a convolved design matrix or a struct of model RDMs.
%
% [res,nulldist,bootdist] = roidata_glm(designvol,epivol,varargin)
function [res,nulldist,bootdist] = roidata_glm(designvol,epivol,varargin)

getArgs(varargin,{'glmclass','GLM','nperm',1,'nboot',0,...
    'glmclassargs',{},'customfun',[],'defaulttail','right',...
    'contrasts',emptystruct('name','conplus','conminus','tail'),...
    'bootmeth','bootstrapruns','bootprep','preparerunboots',...
    'permmeth','permuteruns','permprep','preparerunperms',...
    'customfits',emptystruct('name','funhand','cons','tail')});

if ~iscell(glmclassargs)
    % this extra bit of flexibility helps with cases where you just want to
    % set a single extra parameter without needing to wrap it in curlies
    glmclassargs = {glmclassargs};
end

if ~iscell(customfun)
    if isempty(customfun)
        customfun = {};
    else
        customfun = {customfun};
    end
end
ncustom = numel(customfun);

if isa(designvol,'Volume')
    assert(isequal(epivol.meta.samples.chunks,...
        designvol.meta.samples.chunks),...
        'mismatched epivol/designvol');
    model = vol2glm(designvol,epivol,glmclass,glmclassargs{:});
    res.rows_contrast = designvol.meta.features.labels';
elseif isrdm(designvol)
    model = array2glm(asrdmvec(designvol),epivol.data,...
        ones(size(epivol.data,1),1),glmclass,glmclassargs{:});
    res.rows_contrast = ascol({designvol.name});
end

ncon = model(1).npredictors;

eyecon = eye(model(1).npredictors,class(model(1).data));

% Get the basic estimates (nb we don't just fit directly because in some
% sub-classes the output transform only happens at this level).
res.b = contrast(model,eyecon);
% NB these are OLS p values so probably quite optimistic (you are unlikely
% to have as many df as volumes). Comparison to pperm can be quite
% informative.
res.tail = repmat({defaulttail},[ncon 1]);
paraok = true;
try
    res.sterr = standarderror(model,eyecon);
    res.t = tmap(model,eyecon);
    res.ppara = pmap(model,eyecon,defaulttail);
catch 
    err = lasterror;
    if ~strcmp(err.identifier,'RSA:noParametricStats')
        rethrow(err);
    end
    paraok = false;
end


if isfield(epivol.meta.features,'nfeatures')
    res.nfeatures = epivol.meta.features.nfeatures;
end
res.cols_roi = epivol.meta.features.names;

% custom voodoo
for c = 1:ncustom
    res.custom{c,1} = feval(customfun{c},model);
end

% permutation p values
pind = feval(permprep,model,nperm);
% NB we permute T here to be concordant with the sign flip method in % roidata_rfx. And because we typically want to threshold T maps anyway.
% but for now also b so we can just cat in the slope/contrasts etc (it's
% difficult to calculate t values for these)
[nulldist.b] = feval(permmeth,model,pind,'contrast',[],eyecon);
if paraok
    [nulldist.t] = feval(permmeth,model,pind,'tmap',[],eyecon);
end

% bootstrap estimates
bind = feval(bootprep,model,nboot);
bootdist.b = feval(bootmeth,model,bind,'contrast',[],eyecon);

% store contrasts before non-model parameters get added
regnames = res.rows_contrast;
% deal with contrasts and custom fits
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
        res.b(end+1,:,:) = customfits(c).funhand(...
            res.b(conind,:),customfits(c));
        res.rows_contrast{end+1} = customfits(c).name;
        res.tail{end+1} = customfits(c).tail;
        if nperm > 1
            b = cell(1,size(nulldist.b,3));
            parfor n = 1:size(nulldist.b,3)
                b{n} = customfits(c).funhand(...
                    nulldist.b(conind,:,n),customfits(c));
            end
            nulldist.b(end+1,:,:) = cat(3,b{:});
        end % nperm > 1
        if nboot > 1
            b = cell(1,size(bootdist.b,3));
            for bo = 1:size(bootdist.b,3)
                b{bo} = customfits(c).funhand(...
                    bootdist.b(conind,:,bo),customfits(c));
            end
            bootdist.b(end+1,:,:) = cat(3,b{:});
        end
    end
end

% add in the contrasts
if ~isempty(contrasts)
    if ischar(contrasts)
        contrasts = feval(contrasts,regnames);
    elseif isa(contrasts,'function_handle')
        contrasts = contrasts(regnames);
    end
    assert(isstruct(contrasts),'contrasts must be a struct');
    if ~isfield(contrasts,'convec')
        for c = 1:numel(contrasts)
            contrasts(c).convec = zeros(1,model(1).npredictors);
            posind = ismember(regnames,contrasts(c).conplus);
            % ensure sum to 1
            contrasts(c).convec(posind) = 1/sum(posind);
            nind = ismember(regnames,contrasts(c).conminus);
            contrasts(c).convec(nind) = -1/sum(nind);
        end
    end
    res.rows_contrast = [res.rows_contrast; ascol({contrasts.name})];
    res.tail = [res.tail; ascol({contrasts.tail})];
    conmat = vertcat(contrasts.convec);
    res.b = [res.b; contrast(model,conmat)];
    nulldist.b = vertcat(nulldist.b,feval(permmeth,model,pind,...
        'contrast',[],conmat));
    bootdist.b = vertcat(bootdist.b,feval(bootmeth,model,bind,...
        'contrast',[],conmat));
    if paraok
        res.sterr = [res.sterr; standarderror(model,conmat)];
        res.t = [res.t; tmap(model,conmat)];
        % for tailed tests we need to loop
        for c = 1:numel(contrasts)
            res.ppara = [res.ppara; pmap(model,contrasts(c).convec,...
                contrasts(c).tail)];
        end
        nulldist.t = vertcat(nulldist.t,feval(permmeth,model,pind,...
            'tmap',[],conmat));
    end
end

bootdist.rows_contrast = res.rows_contrast;
bootdist.cols_roi = res.cols_roi;
nulldist.rows_contrast = res.rows_contrast;
nulldist.cols_roi = res.cols_roi;

if nperm > 1
    res.pperm = NaN(size(res.b));
    for ttarget = {'right','left','both'}
        tstr = ttarget{1};
        tind = strcmp(res.tail,tstr);
        res.pperm(tind,:,:) = permpvalue(nulldist.b(tind,:,:,:),tstr);
    end
end
if nboot > 0
    [res.medianboot,res.medianste] = bootprctile(bootdist.b);
end

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
