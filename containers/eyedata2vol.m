% Convert some datastruct to pilab volumes. The output eyevol has one row
% per trial and one feature per sample and modality (any combination of
% x,y,pupil), while the designvol contains dummy variables coding each
% predictor in conind. Used mainly to support wrappers for particular eye
% trackers (see mreyelink2vol and smieye2vol).
%
% INPUT         DEFAULT         DESCRIPTION
% datastruct    -               struct with fields x,y,pupil (each ntrial
%                                   by nsample), time (1 by nsample),conind
%                                   (ntrial by 1)
% frameperiod   -               sample time in s
%
% NAMED INPUT
% targetfields  {'x','y','pupil'}   select data channels for import.
% demean        false           de-mean each run
% dur           -               duration of trials in s. Defaults to
%                                   keeping all samples in each
%                                   targetfield.
% downperiod    []              downsample signal to this rate
% removetrialthresh     .5      remove trials with more than this much
%                                   missing data
% screenlims    -               rescore data outside these lims as missing
%
% [eyevol,designvol] = eyedata2vol(datastruct,frameperiod,varargin)
function [eyevol,designvol] = eyedata2vol(datastruct,frameperiod,varargin)

getArgs(varargin,{'demean',false,'downperiod',[],...
    'removetrialthresh',.5,'screenlims',[],'dur',[],...
    'targetfields',{'x','y','pupil'}});


dodownsample = ~isempty(downperiod) && downperiod~=frameperiod;
[ntrials,npossiblesamples] = size(datastruct.x);
% if you don't specify a duration we assume you want it all
nsamples = npossiblesamples;
if ~isempty(dur)
    nsamples = round(dur/frameperiod);
end
% keep track of bad trials (missing data or ignorelabels)
badtrials = false([ntrials,1]);
assert(nsamples <= npossiblesamples,'not enough samples for dur');

% build design matrix
assert(ndims(datastruct.conind)==2,...
    'conind must be 2D design matrix or vector of indices')
if all(size(datastruct.conind)>1)
    % assume pre-cooked design matrix
    design = datastruct.conind;
else
    ulabels = uniquen(datastruct.conind);
    linind = sub2ind([ntrials numel(ulabels)],ascol(1:ntrials),...
        ascol(datastruct.conind));
    design = zeros(ntrials,numel(ulabels));
    design(linind) = 1;
end

% rescore data outside screen edges as missing
if ~isempty(screenlims)
    outsideind = abs(datastruct.x) > screenlims(1) | ...
        abs(datastruct.y)>screenlims(2);
    datastruct.x(outsideind) = NaN;
    datastruct.y(outsideind) = NaN;
    datastruct.pupil(outsideind) = NaN;
end

% preprocess data, figure out which trials are bad
ndown = nsamples;
for target = targetfields(:)'
    targetstr = target{1};
    % restrict to intended number of samples
    datastruct.(targetstr) = datastruct.(targetstr)(:,1:nsamples);

    if dodownsample
        % downsample the signal to this many samples
        ndown = round(nsamples * (frameperiod / downperiod));
        % raw indices for each new sample
        downind = floor(linspace(1,ndown+1,nsamples+1))';
        downind(end) = [];
        % so now we can accumarray based on this
        newdata = NaN([ntrials ndown],class(datastruct.(targetstr)));
        % nan to interpolate missing data whenever possible, median to
        % minimise outlier influence
        for tr = 1:ntrials
            newdata(tr,:) = accumarray(downind,...
                datastruct.(targetstr)(tr,:),[],@nanmedian);
        end
        % put back in datastruct
        datastruct.(targetstr) = newdata;
    end

    % update badtrials with trials that have too much missing data
    % (suggesting interpolation can't be trusted)
    % (but note that we keep it in for now because it's easier to just
    % do all the bad data removal at once at the end)
    pmissing = sum(isnan(datastruct.(targetstr)),2) ./ ...
        size(datastruct.(targetstr),2);
    badtrials = badtrials | pmissing > removetrialthresh;

    % do the interpolation (transpose since interpolatemissing works on
    % rows not columns)
    datastruct.(targetstr) = interpolatemissing(...
        datastruct.(targetstr)')';

    if demean
        % now we need to de-mean relative to the INTERPOLATED mean. Bit
        % funky.
        datastruct.(targetstr) = datastruct.(targetstr) - nanmean(...
            datastruct.(targetstr)(:));
    end
end
if dodownsample
    timeind = [true; diff(downind)~=0];
    datastruct.time = datastruct.time(timeind);
end


eyevol = {};
designvol = {};
if all(badtrials)
    logstr('no valid data remains, returning empty\n')
    return
end

% second pass to remove bad trials and construct eyevol
design(badtrials,:) = [];
for target = targetfields(:)'
    targetstr = target{1};
    datastruct.(targetstr)(badtrials,:) = [];
    thismeta = struct('time',asrow(datastruct.time),'datatype',...
        {repmat({targetstr},[1 ndown])});
    eyevol{end+1} = Volume(datastruct.(targetstr),'metafeatures',...
        thismeta);
end
eyevol = horzcat(eyevol{:});
designvol = Volume(design,'metafeatures',struct('labels',...
        {mat2strcell(1:size(design,2),['conind_%d'])}));
