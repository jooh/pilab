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
% scans         1:viewGet(v,'nscans')   scans to import
% ignorelabels  []              remove data for these values of varname
% taskNum       1               inputs for getTaskEyeTraces
% phaseNum      1
% segmentNum    1
% dataPad       3
% removeBlink   .15             extra time to cut around each blink to
%                                   reduce eyelid artefact
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
ulabels = [];

rc = 0;
eyevol = {};
designvol = {};
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
    % check that labels aren't changing over runs
    ulabels = setifunset(ulabels,unique(vardata),true);

    % TODO: also get conind
    datastruct = struct('x',t.eye.xPos,'y',t.eye.yPos,...
        'pupil',t.eye.pupil,'time',t.eye.time,'conind',vardata);

    % keep track of how much we're missing per run
    pmissing = sum(asrow(isnan(datastruct.x(:,1:nsamples)))) ./ ...
        numel(datastruct.x(:,1:nsamples));
    if pmissing > removerunthresh
        logstr('excessive missing data, skipping.\n';
        continue
    end
    % so we're keeping this run
    rc = rc+1;
    [eyevol{rc},designvol{rc}] = eyedata2vol(datastruct,frameperiod,...
        'screenlims',screenlims,...
        'downperiod',downperiod,...
        'demean',demean,...
        'removetrialthresh',removetrialthresh,...
        'dur',dur,...
        'targetfields',targetfields);
    ntrials = eyevol{rc}.nsamples;
    eyevol{rc}.meta.samples.chunks = ones(ntrials,1)*r;
    designvol{rc}.meta.samples.chunks = ones(ntrials,1)*r;
    eyevol{rc}.meta.samples.filenames = repmat(stimfiles(r),[ntrials 1]);
    designvol{rc}.meta.samples.filenames = repmat(stimfiles(r),[ntrials 1]);
end
