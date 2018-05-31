% apply standard filtering and preprocessing options to an datavol/designvol
% pair.
%
% INPUT             DEFAULT     DESCRIPTION
% datavol            -          Volume instance
% designvol         -           Volume instance
%
% NAMED INPUT       DEFAULT     DESCRIPTION
% matchn            false       match number of volumes across runs
% covariatedeg      ''          degree of polynomial detrend (0:n).
%                                   'adaptive' means we find a value using
%                                   vol2covdeg
% domedianfilter    false       medianfilter the data over time
% sgdetrend         false       Savitzky-Golay detrend the data over time
% sgolayK           NaN         parameter for Savitzky-Golay detrend
% sgolayF           NaN         parameter for Savitzky-Golay detrend
% dozscore          false       zscore the data over time
% targetlabels      []          cell array containing names for any labels
%                                   in designvol to select (keep all if
%                                   isempty)
% ignorelabels      []          cell array containing names for any labels
%                                   to remove
% setclass          []          convert all data to specified class
% resortind         []          resort the order of designvol according to
%                                   these numerical indices
% percentsignal     false       if true, convert data to % signal change
%                                   and scale peak of design matrix to 1.
%                                   Note that this doesn't perfectly work
%                                   for rapid ER designs.
%
% [datavol,designvol] = preprocessvols(datavol,designvol,varargin)
% Modified by TEC 17/5/18 to allow multiple ignoreinds as a cell array
% Modified by JC and TEC 31/5/18 to allow design volume changes only

function [datavol,designvol] = preprocessvols(datain,designin,varargin)

getArgs(varargin,{'matchn',false,'covariatedeg','',...
    'domedianfilter',false,'sgdetrend',false,'sgolayK',NaN,'sgolayF',NaN,...
    'dozscore',false,'targetlabels',[],'ignorelabels',[],...
    'setclass',[],'resortind',[],'percentsignal',false});

% need to make a copy of each volume here, otherwise changes are made in
% place, quietly, in some cases.

hasdata = true;
if ~exist('datain','var') || isempty(datain)
    hasdata = false;
end

if hasdata
    datavol = datain(:,1:end);
    assert(strcmp(class(datavol),class(datain)),'indexing broke data class')
end
designvol = designin(:,1:end);
% this test shouldn't be here but there have historically been problems
% with this kind of stuff so let's just be extra sure.

assert(strcmp(class(designvol),class(designin)),...
    'indexing broke data class')

if matchn
    nperchunk = arrayfun(@(c)sum(datavol.meta.samples.chunks==c),...
        datavol.desc.samples.unique.chunks);
    targetn = min(nperchunk);
    if all(targetn == nperchunk);
        logstr('all chunks have same number of samples\n');
    else
        logstr('matching nsamples to smallest chunk (%d)\n',...
            targetn);
        % this might actually be a tad hairy. need to apply the
        % same to design and epi obviously
        goodsamp = false(datavol.nsamples,1);
        for c = 1:datavol.desc.samples.nunique.chunks
            chunkind = find(datavol.meta.samples.chunks == ...
                datavol.desc.samples.unique.chunks(c));
            goodsamp(chunkind(1:targetn)) = true;
        end
        datavol = datavol(goodsamp,:);
        designvol = designvol(goodsamp,:);
        logstr('removed %d samples (%2.0f%% of total)\n',...
            sum(~goodsamp),100*sum(~goodsamp) / numel(goodsamp));
    end
end

if percentsignal
    logstr('scaling data to percent signal change\n'); 
    if hasdata
        %scale epi by 100 * mean per voxel
        filterbychunk(datavol,@(x)100*bsxfun(@rdivide,x,mean(x)));
    end
    % scale design matrix by global max (so the peak across the design
    % matrix is 1)
    filterbychunk(designvol,@(x)bsxfun(@rdivide,x,max(x(:))));
end

% de-trend config
if strcmp(covariatedeg,'adaptive')
    % does this work without datavol?
    covariatedeg = vol2covdeg(designvol);
end

% first high-pass trend removal
if ~isempty(covariatedeg)
    logstr('polynomial detrend (degree=%.0f)\n',covariatedeg);
    if hasdata
        filterbychunk(datavol,'polydetrend',covariatedeg);
    end
    filterbychunk(designvol,'polydetrend',covariatedeg);
end

if ~isempty(domedianfilter) && domedianfilter
    logstr('median filter (n=%.0f)\n',domedianfilter);
    if hasdata
        filterbychunk(datavol,'medianfilter',domedianfilter);
    end
    filterbychunk(designvol,'medianfilter',domedianfilter);
end

if sgdetrend
    logstr('Savitzky-Golay detrend (k=%.0f,f=%.0f)\n',...
        sgolayK,sgolayF);
    % ensure double
    if hasdata
        datavol.data = double(datavol.data);
        sgdetrend(datavol,sgolayK,sgolayF);
    end
    designvol.data = double(designvol.data);    
    sgdetrend(designvol,sgolayK,sgolayF);
end

if dozscore
    logstr('Z-scoring samples\n')
    if hasdata
        filterbychunk(datavol,'zscore',[],1);
    end
    filterbychunk(designvol,'zscore',[],1);
end

% handle selecting labels
if ~isempty(targetlabels)
    [~,~,coninds] = intersect(targetlabels, designvol.meta.features.labels, ...
        'stable');
    assert(~isempty(coninds),...
    'found no labels matching targetlabels %s',cell2str(targetlabels));
    designvol = designvol(:,coninds);
end

% handle removing and projecting out labels
if ~isempty(ignorelabels)
    % things to project out
    [~,covinds] = intersect(designvol.meta.features.labels, ...
        ignorelabels,'stable');
    assert(~isempty(covinds), 'no labels matching ignorelabels %s', ...
        cell2str(ignorelabels));
    
    % things to keep
    [~,keepinds] = setdiff(designvol.meta.features.labels, ...
        ignorelabels,'stable');
    
    covariates = designvol.data(:,covinds);
    designvol = designvol(:,keepinds);
    logstr('projecting out %d covariates\n',size(covariates,2));
    for c = 1:designvol.desc.samples.nunique.chunks
        chunkind = designvol.meta.samples.chunks == ...
            designvol.desc.samples.unique.chunks(c);
        if hasdata
            assert(isequal(chunkind,...
                datavol.meta.samples.chunks == ...
                designvol.desc.samples.unique.chunks(c)),...
                'mismatched chunks in datavol and designvol');
            datavol.data(chunkind,:) = projectout(...
                datavol.data(chunkind,:),covariates(chunkind,:));
        end
        designvol.data(chunkind,:) = projectout(...
            designvol.data(chunkind,:),covariates(chunkind,:));
    end
end

% now set class last - so maximal precision for pre-processing
if ~isempty(setclass)
    logstr('setting data to %s\n',setclass);
    if hasdata
        datavol.data = feval(setclass,datavol.data);
    end
    designvol.data = feval(setclass,designvol.data);
end

if ~isempty(resortind)
    if ischar(resortind)
        resortind = feval(resortind,designvol.nfeatures);
    end
    logstr('resorting regressors in designvol.\n');
    designvol = designvol(:,resortind);
end

if ~hasdata
    datavol = [];
end
