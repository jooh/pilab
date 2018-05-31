% Convert samples txt file from SMI eye tracker to pilab volume instances.
%
% INPUT         DEFAULT         DESCRIPTION
% filenames     -               cell array of paths to text files from SMI
%                                   samples converter
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
% ignorelabels  []              remove data for these values of varname
% removeBlink   .15             extra time to cut around each blink to
%                                   reduce eyelid artefact
%
% [eyevol,designvol] = smieye2vol(stimfiles,varname,dur,varargin)
%
function [eyevol,designvol] = smieye2vol(filenames,dur,varargin)

getArgs(varargin,{'targetfields',{'x','y','pupil'},'demean',false,...
    'downperiod',[],'removerunthresh',.25,'removetrialthresh',.5,...
    'ignorelabels',[],'removeBlink',.3});

if ~iscell(filenames)
    filenames = {filenames};
end

nfile = numel(filenames);
expar = struct;
rc = 0;
for f = 1:nfile
    [samples,expar] = em_importsamples(filenames{f},expar);
    % convert ms to s
    samples(:,1) = samples(:,1) / 1000;
    % now convert to eyemat
    % so in this awkward format, samples(:,2) contains numerical codes for
    % trial labels (see expar.markers), samples(:,3) is pupil diameter,
    % (:,4) is X and (:,5) is Y. 
    assert(issorted(samples(:,1)))
    % work out frameperiod empirically
    rawdata = samples(isnan(samples(:,2)),[1 3:5]);
    sr = diff(rawdata(:,1));
    frameperiod = median(diff(rawdata(:,1)));
    nsamples = round(dur/frameperiod);
    % make sure any zeros (blinks) are rescored to nan
    blinks = any(rawdata(:,2:end)==0,2);
    % so here we need to be careful to index only cols 2:4
    rawdata([false(size(rawdata,1),1) repmat(blinks,[1 3])]) = NaN;
    % pad blinks
    if removeBlink
        % remove period in sample units
        nrem = round(removeBlink/frameperiod);
        % find current blinks
        missing = single(isnan(rawdata(:,2:end)));
        % convolution is a convenient way to expand every missing sample to
        % nrem.
        smcheck = conv2(missing,repmat(1/nrem,nrem,1),'same');
        rawdata([false([size(rawdata,1),1]),smcheck > 0]) = NaN;
    end
    messages = samples(~isnan(samples(:,2)),1:2);
    % remove ignored labels (e.g. SOA)
    for lab = ignorelabels(:)'
        labstr = lab{1};
        messages(messages(:,2) == expar.markers.(labstr),:) = [];
    end
    % rescore to numerical index in alpha order
    labels = sort(fieldnames(expar.markers));
    labels(strcmp(labels,'Fixation')) = [];
    labels(strcmp(labels,'Saccade')) = [];
    labels(strcmp(labels,'Blink')) = [];
    for lab = 1:numel(labels)
        mcode = expar.markers.(labels{lab});
        messages(messages(:,2)==mcode,2) = lab;
    end
    ntrials = size(messages,1);
    data.x = NaN([ntrials,nsamples]);
    data.y = NaN([ntrials,nsamples]);
    data.pupil = NaN([ntrials,nsamples]);
    data.conind = messages(:,2);
    data.time = NaN([ntrials,nsamples]);
    for t = 1:ntrials
        starttime = messages(t,1);
        endtime = starttime + dur;
        % first data sample after starttime
        firstind = find(rawdata(:,1)>starttime,1,'first');
        lastind = find(rawdata(:,1)<=endtime,1,'last');
        trialdata = rawdata(firstind:lastind,2:4);
        ntrialsamples = size(trialdata,1);
        data.pupil(t,1:ntrialsamples) = trialdata(:,1);
        data.x(t,1:ntrialsamples) = trialdata(:,2);
        data.y(t,1:ntrialsamples) = trialdata(:,3);
        data.time(t,1:ntrialsamples) = (rawdata(firstind:lastind,1) - ...
            rawdata(firstind,1))';
    end
    % reduce precision of time dim to avoid spurious mismatch errors when
    % concatenating runs.
    data.time = reduceprecision(mean(data.time),3);
    pmissing = sum(asrow(isnan(data.x))) ./ numel(data.x);
    if pmissing > removerunthresh
        logstr('excessive missing data, skipping.\n');
        continue
    end
    rc = rc+1;
    [eyevol{rc},designvol{rc}] = eyedata2vol(data,frameperiod,...
        'demean',demean,'downperiod',downperiod,...
        'removetrialthresh',removetrialthresh,...
        'screenlims',expar.screen_res_px,'dur',dur,...
        'targetfields',targetfields);
    if isempty(eyevol{rc})
        rc = rc-1;
        logstr('excessive missing trials, skipping.\n');
        continue
    end

    ntrials = eyevol{rc}.nsamples;
    eyevol{rc}.meta.samples.chunks = ones(ntrials,1)*f;
    designvol{rc}.meta.samples.chunks = ones(ntrials,1)*f;
    eyevol{rc}.meta.samples.filenames = repmat(filenames(f),[ntrials 1]);
    designvol{rc}.meta.samples.filenames = repmat(filenames(f),[ntrials 1]);
end
% finally, concatenate the whole thing
eyevol = cat(1,eyevol{:});
designvol = cat(1,designvol{:});

% Import samples for data analysis. Optionally enter a parameters file to
% produce consistent trigger coding across datasets - also checks that
% other parameters are consistent across datasets.
% [samples,parameters] = em_parsesamplestxt(samplests,[parameters])
function [samples, par] = em_importsamples(samplestxt,expar)

% Parse the sample text file
% [samples, par] = em_parsesamplestxt(samplestxt)
fid = fopen(samplestxt);
% Read header info
nhlines = 47;
theader = textscan(fid,'%s',nhlines,'delimiter','\n','whitespace','');
theader = theader{1};
% Read data
tdata = textscan(fid,'%s','delimiter','\n','whitespace','');
tdata = tdata{1};
fclose(fid);
% Now to make sense of things...

% Extract parameters from this session
datetime = textscan(theader{3}, '%*s %*s %s %s');
par.rec_date = datetime{1}{1};
par.rec_time = datetime{2}{1};
txtpath = textscan(theader{2},'%*s %*s %*s %s');
par.txt_path = txtpath{1}{1};
samplingr = textscan(theader{6},'%*s %*s %*s %f');
par.sampling_rate_hz = samplingr{1};
filename = textscan(theader{13},'%*s %*s %s');
par.filename = filename{1}{1};
screenres = textscan(theader{17},'%*s %*s %*s %f %f');
par.screen_res_px = cell2mat(screenres);
screenwidth = textscan(theader{32},'%*s %*s %*s %*s %f %*d');
par.screenwidth = screenwidth{1};
viewdistance = textscan(theader{33},'%*s %*s %*s %*s %f');
par.viewdistance = viewdistance{1};
% Detect binocular acquisitions
par.binocular = ~isempty(strfind(theader{45},'RIGHT'));

% Marker code start number
mc = 10e3;

% No previous markers to consider, of checks to make
if ~exist('expar','var') || ~isfield(expar,'markers')
	% Initialise empty markers struct
	expar.markers = struct;
end

% Make some basic consistency checks
for fn = {'sampling_rate_hz','screenwidth','viewdistance','binocular'}
	fn = fn{1};
	if isfield(expar,fn)
		if ~par.(fn) == expar.(fn)
			error('Inconsistent parameters for %s %s', ...
				par.filename,expar.filename)
		end
	end
end

% Figure out how many markers we've already used
codes = cell2mat(structfun(@(x){x},expar.markers));
% And ensure we don't duplicate
mc = max([max(codes) mc]);

% Now update par with expar
par = catstruct(par,expar);

% NB These fields are not actually used now - just reserved for later use
% in event extraction
if ~isfield(par.markers,'Saccade');
	mc = mc + 1e3;
	par.markers.Saccade = mc;
	mc = mc + 1e3;
	par.markers.Blink = mc;
	mc = mc + 1e3;
	par.markers.Fixation = mc;
end

% Initialise output matrix
nsamples = length(tdata);
% row: [timestamp eventflag Lx Ly [Rx Ry]]

% textscan setup for messages
msgstr = '%f %*s %*d %*s %*s %s';
if par.binocular
    error('currently unsupported')
else
	samples = zeros(nsamples,4);
	% textscan for monocular data
	scanstr = '%f %*s %*d %*f %*f %f %*f %*f %*f %f %f %*f';
	linewidth = 5;
end
samples(:) = NaN;

% Now parse samples
for s = 1:length(tdata);
	% Behaviour depends on line type
	switch isempty(strfind(tdata{s},'MSG'));
		case 1 % regular data line (SMP)
			% 3 or 5 entries
			lineout = cell2mat(textscan(tdata{s},scanstr));
			% Column indices to insert data into
			lineIs = [1 3:linewidth];
		case 0 % message line
			linedata = textscan(tdata{s},msgstr);
			% Strip any extensions
			[x,fn,x] = fileparts(linedata{2}{1});
            if ~isfield(par.markers,fn)
                mc = mc + 1e3;
                par.markers.(fn) = mc;
            end
            lineout = [linedata{1} par.markers.(fn)];
			lineIs = [1 2];
	end
	% Add data to outmat
	samples(s,lineIs) = lineout;
end

% Detect sorting errors in sample log
if ~issorted(samples(:,1))
	fprintf('Repaired sample sort error in %s', ...
		par.filename)
	[x,I] = sort(samples(:,1));
	samples = samples(I,:);
end

% Rescore time to ms
samples(:,1) = samples(:,1) * .001;
