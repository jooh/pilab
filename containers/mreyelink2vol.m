% Convert mrTools / MGL eyelink eye tracker data to pilab volume instances.
% The data for each trial is concatenated horizontally to make a feature
% vector. See also mrtools2vol.
%
% INPUT         DEFAULT         DESCRIPTION
% stimfiles     -               cell array of stim mats or a view handle
% varname       -               variable name to target in stimfiles
% dur           -               duration of trials in s
%
% NAMED INPUT
% targetfields  {'x','y','pupil'}   select data channels for import.
% demean        false           de-mean each run
% downperiod    []              downsample signal to this rate
% removerunthresh   .25         remove runs with more than this much
%                                   missing data
% removetrialthresh     .5      remove trials with more than this much
%                                   missing data
% screenlims    d.myscreen.displaySize/2 rescore data outside these lims as
%                                   missing
%
% scans         1:viewGet(v,'nscans')   scans to import
% ignorelabels  []              remove data for these values of varname
% taskNum       1               inputs for getTaskEyeTraces
% phaseNum      1
% segmentNum    1
% dataPad       3
% removeBlink   .15
%
% [eyevol,designvol] = mreyelink2vol(stimfiles,varname,dur,varargin)
function [eyevol,designvol] = mreyelink2vol(stimfiles,varname,dur,varargin)

getArgs(varargin,{'taskNum',1,'phaseNum',1,'segmentNum',1,...
    'removeBlink',.15,'dataPad',3,'targetfields',[],'demean',false,...
    'downperiod',[],'removerunthresh',.25,'removetrialthresh',.5,...
    'screenlims',[],'scans',[],'ignorelabels',[]});

% parse input
if isview(stimfiles)
    v = stimfiles;
    logstr('extracting stimfiles from %s\n',viewGet(v,'groupname'));
    clear stimfiles
    if isempty(scans)
        scans = 1:viewGet(v,'nscans');
    end
    for r = 1:numel(scans)
        thisscan = scans(r);
        s = viewGet(v,'stimfilename',scans(r));
        assert(numel(s)==1,...
            'each scan must be associated with a single stim file');
        stimfiles(r) = s;
    end
end
nrun = numel(stimfiles);

% map from mrtools to own variable names
fieldmap = cell2struct({'xPos','yPos','pupil'},{'x','y','pupil'},2);

% possibly use only some fields
if ~isempty(targetfields)
    badfield = setdiff(fieldnames(fieldmap),targetfields);
    fieldmap = rmfield(fieldmap,badfield);
end

frameperiod = [];

metasamples = struct('chunks',[],'filenames',[]);
datametafeatures = struct('time',[]);
ulabels = [];

rc = 0;
for r = 1:nrun
    t = getTaskEyeTraces(stimfiles{r},'dataPad',dataPad,'removeBlink',...
        removeBlink,'taskNum',taskNum,'phaseNum',phaseNum,'segNum',...
        segmentNum);
    if ~isfield(t,'eye')
        logstr('no eye data, skipping run %d\n',r);
        continue;
    end

    if isempty(screenlims)
        % infer from myscreen
        d = load(stimfiles{r});
        screenlims = d.myscreen.displaySize/2;
    end

    frameperiod = setifunset(frameperiod,median(diff(t.eye.time)),1);
    dodownsample = ~isempty(downperiod) && downperiod~=frameperiod;
    nsamples = round(dur/frameperiod);
    [ntrials,npossiblesamples] = size(t.eye.xPos);
    % keep track of bad trials (missing data or ignorelabels)
    badrows = false([ntrials,1]);
    assert(nsamples <= npossiblesamples,'not enough samples for dur');

    % rescore data outside screen edges as missing
    outsideind = abs(t.eye.xPos) > screenlims(1) | ...
        abs(t.eye.yPos)>screenlims(2);
    t.eye.xPos(outsideind) = NaN;
    t.eye.yPos(outsideind) = NaN;
    t.eye.pupil(outsideind) = NaN;

    % keep track of how much we're missing per run
    pmissing = sum(asrow(isnan(t.eye.xPos(:,1:nsamples)))) ./ ...
        numel(t.eye.xPos(:,1:nsamples));
    if pmissing > removerunthresh
        logstr('excessive missing data, skipping run %d\n',r);
        continue;
    end
    rc = rc+1;
    data.qa(rc).propmissing = pmissing;
    assert(isequal(isnan(t.eye.xPos),isnan(t.eye.yPos),...
        isnan(t.eye.pupil)),'mismatched NaN across data fields');
    data.qa(rc).wasmissing = isnan(t.eye.xPos);

    % collect data
    thistime = t.eye.time(1:nsamples);
    for outfn = fieldnames(fieldmap)'
        outfnstr = outfn{1};
        mrfnstr = fieldmap.(outfnstr);
        data.raw(rc).(outfnstr) = t.eye.(mrfnstr)(:,1:nsamples);

        if dodownsample
            % downsample the signal to this many samples
            ndown = round(nsamples * (frameperiod / downperiod));
            % raw indices for each new sample
            downind = floor(linspace(1,ndown+1,nsamples+1))';
            downind(end) = [];
            % so now we can accumarray based on this
            newdata = NaN([ntrials ndown],class(data.raw(rc).(outfnstr)));
            for tr = 1:ntrials
                newdata(tr,:) = accumarray(downind,...
                    data.raw(rc).(outfnstr)(tr,:),[],@nanmedian);
            end
            % put back in data.raw
            data.raw(rc).(outfnstr) = newdata;
        end

        % extract means before interpolation, but after downsampling
        data.means(rc).(outfnstr) = nanmean(data.raw(rc).(outfnstr),2);
        
        % update badrows with trials that have too much missing data
        % (suggesting interpolation can't be trusted)
        % (but note that we keep it in for now because it's easier to just
        % do all the bad data removal at once at the end)
        pmissing = sum(isnan(data.raw(rc).(outfnstr)),2) ./ ...
            size(data.raw(rc).(outfnstr),2);
        badrows = badrows | pmissing > removetrialthresh;

        % do the interpolation (transpose since interpolatemissing works on
        % rows not columns)
        data.raw(rc).(outfnstr) = interpolatemissing(...
            data.raw(rc).(outfnstr)')';

        if demean
            % now we need to de-mean relative to the INTERPOLATED mean. Bit
            % funky.
            data.raw(rc).(outfnstr) = data.raw(rc).(outfnstr) - mean(...
                data.raw(rc).(outfnstr)(:));
        end
    end
    if dodownsample
        timeind = [true; diff(downind)~=0];
        thistime = thistime(timeind);
    end

    % extract meta data
    metasamples(rc).chunks = ones(ntrials,1)*r;
    metasamples(rc).filenames = repmat(stimfiles(r),[ntrials 1]);
    datametafeatures.time = setifunset(datametafeatures.time,thistime);

    % figure out model
    if isfield(t.parameter,varname)
        % use that
        d = t.parameter;
    elseif isfield(t.randVars,varname)
        % use that
        d = t.randVars;
    else
        error('could not find varname: %s',varname);
    end
    vardata = d.(varname);
    if ~isempty(ignorelabels)
        badrows = badrows | vardata' == ignorelabels;
        vardata(badrows) = [];
        for fn = fieldnames(data.raw)'
            fnstr = fn{1};
            data.raw(rc).(fnstr)(badrows,:) = [];
        end
        for fn = fieldnames(metasamples)'
            fnstr = fn{1};
            metasamples(rc).(fnstr)(badrows) = [];
        end
        % update ntrials estimate
        ntrials = length(metasamples(rc).(fnstr));
    end
    ulabels = setifunset(ulabels,unique(vardata),true);
    linind = sub2ind([ntrials numel(ulabels)],1:ntrials,vardata);
    design{rc} = zeros(ntrials,numel(ulabels));
    design{rc}(linind) = 1;
    if demean
        design{rc} = bsxfun(@minus,design{rc},mean(design{rc},1));
    end
end

designvol = [];
eyevol = {};
if rc>0
    % NB frameperiod isn't strictly correct here - there may be gaps between
    % trials
    metasamples = collapsestruct(metasamples,@vertcat,true);
    designvol = Volume(cat(1,design{:}),'metasamples',metasamples,...
        'metafeatures',struct('labels',...
        {mat2strcell(ulabels,[varname '_%d'])}));

    % now combine some set of the data and make the eyevol
    for outfn = fieldnames(fieldmap)'
        outfnstr = outfn{1};
        thismeta = datametafeatures;
        thismeta.datatype = repmat({outfnstr},[1 size(data.raw(1).(outfnstr),2)]);
        eyevol{end+1} = Volume(cat(1,data.raw.(outfnstr)),'metasamples',...
            metasamples,'metafeatures',thismeta);
    end
    % and collapse to a single eyevol
    eyevol = horzcat(eyevol{:});
end
