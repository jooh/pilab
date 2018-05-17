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
function [datavol,designvol] = preprocessvols(datain,designin,varargin)

getArgs(varargin,{'matchn',false,'covariatedeg','',...
    'domedianfilter',false,'sgdetrend',false,'sgolayK',NaN,'sgolayF',NaN,...
    'dozscore',false,'targetlabels',[],'ignorelabels',[],...
    'setclass',[],'resortind',[],'percentsignal',false});

% need to make a copy of each volume here, otherwise changes are made in
% place, quietly, in some cases.
datavol = datain(:,1:end);
designvol = designin(:,1:end);
% this test shouldn't be here but there have historically been problems
% with this kind of stuff so let's just be extra sure.
assert(strcmp(class(datavol),class(datain)),'indexing broke data class')
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
    % scale epi by 100 * mean per voxel
    filterbychunk(datavol,@(x)100*bsxfun(@rdivide,x,mean(x)));
    % scale design matrix by global max (so the peak across the design
    % matrix is 1)
    filterbychunk(designvol,@(x)bsxfun(@rdivide,x,max(x(:))));
end

% de-trend config
if strcmp(covariatedeg,'adaptive')
    covariatedeg = vol2covdeg(datavol);
end

% first high-pass trend removal
if ~isempty(covariatedeg)
    logstr('polynomial detrend (degree=%.0f)\n',covariatedeg);
    filterbychunk(datavol,'polydetrend',covariatedeg);
    filterbychunk(designvol,'polydetrend',covariatedeg);
end

if ~isempty(domedianfilter) && domedianfilter
    logstr('median filter (n=%.0f)\n',domedianfilter);
    filterbychunk(datavol,'medianfilter',domedianfilter);
    filterbychunk(designvol,'medianfilter',domedianfilter);
end

if sgdetrend
    logstr('Savitzky-Golay detrend (k=%.0f,f=%.0f)\n',...
        sgolayK,sgolayF);
    % insure double
    datavol.data = double(datavol.data);
    designvol.data = double(designvol.data);
    sgdetrend(datavol,sgolayK,sgolayF);
    sgdetrend(designvol,sgolayK,sgolayF);
end

if dozscore
    logstr('Z-scoring samples\n')
    filterbychunk(datavol,'zscore',[],1);
    filterbychunk(designvol,'zscore',[],1);
end

% maybe deal with covariates
nlabels = length(designvol.desc.features.unique.labels);
if isempty(targetlabels)
    coninds = 1:nlabels;
else
    coninds = find(strcmp(designvol.desc.features.unique.labels,...
        targetlabels));
end
if ~isempty(ignorelabels)
    if iscell
        for this_label = 1:length(ignorelabels)
            ignoreinds(this_label) = find(strcmp(designvol.meta.features.labels,ignorelabels{this_label}));
        end
    else
    ignoreinds = find(strcmp(designvol.meta.features.labels,...
        ignorelabels));
    end
    coninds = setdiff(coninds,ignoreinds,'stable');
end
assert(~isempty(coninds),...
    'found no labels matching targetlabels/ignorelabels %s/%s',...
    targetlabels,ignorelabels);
if ~isequal(coninds,1:nlabels)
    covariates = designvol.data(:,ignoreinds);
    designvol = designvol(:,coninds);
    % project out bad cons
    if ~isempty(ignoreinds)
        logstr('projecting out %d covariates\n',...
            size(covariates,2));
        for c = 1:datavol.desc.samples.nunique.chunks
            chunkind = datavol.meta.samples.chunks == ...
                datavol.desc.samples.unique.chunks(c);
            assert(isequal(chunkind,...
                designvol.meta.samples.chunks == ...
                designvol.desc.samples.unique.chunks(c)),...
                'mismatched chunks in datavol and designvol');
            datavol.data(chunkind,:) = projectout(...
                datavol.data(chunkind,:),covariates(chunkind,:));
            designvol.data(chunkind,:) = projectout(...
                designvol.data(chunkind,:),covariates(chunkind,:));
        end
    end
end

% now set class last - so maximal precision for pre-processing
if ~isempty(setclass)
    logstr('setting data to %s\n',setclass);
    datavol.data = feval(setclass,datavol.data);
    designvol.data = feval(setclass,designvol.data);
end

if ~isempty(resortind)
    if ischar(resortind)
        resortind = feval(resortind,designvol.nfeatures);
    end
    logstr('resorting regressors in designvol.\n');
    designvol = designvol(:,resortind);
end
