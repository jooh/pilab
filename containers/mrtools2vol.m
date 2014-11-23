% convert data for a mrTools view to volume instances. Data are loaded in
% 2d masked form. This is basically an intermediate data stage before
% further conversion to GLM instance for decoding.
%
% INPUT         DEFAULT         DESCRIPTION
% v             newView         mrView session
% mask                          the name of an ROI in the subject's ROIs
%                                   dir (e.g. from tseries_roimask)
% varname                       input for getStimvol
% duration                      duration of events in volume units
%
% NAMED INPUT
% group         'MotionComp'    group to target 
% scans         'all'           scan numbers (can be vector instead)
% ignorelabels  []              any values of varname to ignore
% taskNum       1               inputs for getStimvol
% phaseNum      1
% segmentNum    1
% hrffun        []              input for convolveonsets
% checkframe    false           if true, infer from stimfile
%
% OUTPUT
% epivol    MriVolume instance containing detrended in-mask EPI time
%   courses
% designvol Volume instance containing the convolved and detrended
%   design matrix
%
% [epivol,designvol] = mrtools2vol([v],mask,varname,dur,[varargin])
function [epivol,designvol] = mrtools2vol(v,mask,varname,dur,varargin)

getArgs(varargin,{'group','MotionComp','scans','all','taskNum',1,...
    'phaseNum',1,'segmentNum',1,'ignorelabels',[],'hrffun',[],...
    'checkframe',false});

madeview = false;
if ieNotDefined('v')
    v = getMLRView;
    if isempty(v)
        v = newView;
        madeview = true;
    end
end

% easy peasy
epivol = MrToolsVolume(v,mask,'group',group,'scans',scans);

if checkframe
    for s = 1:epivol.desc.samples.nunique.chunks
        thisscan = epivol.desc.samples.unique.chunks(s);
        stimfile = viewGet(v,'stimfile',thisscan);
        volEvents = find(stimfile{1}.myscreen.events.tracenum == 1);
        volTrigRatio = viewGet(v,'auxParam','volTrigRatio',thisscan);
        frameperiod(s) = median(diff(...
            stimfile{1}.myscreen.events.time(volEvents))) / volTrigRatio;
    end
    mframe = mean(frameperiod);
    if mframe ~= epivol.frameperiod
        logstr('estimated frameperiod=%.2f(+/-%.2f) (was %.2f)\n',...
            mframe,std(frameperiod),epivol.frameperiod);
    end
    epivol.frameperiod = mframe;
end

% now obtain the design information in the format required by
% convolveonsets
% first need to re-open the view since MrToolsVolume selfishly destroys it
v = newView;
v = viewSet(v,'currentgroup',group);
X = {};
for s = 1:epivol.desc.samples.nunique.chunks
    thisscan = epivol.desc.samples.unique.chunks(s);
    v = viewSet(v,'currentscan',thisscan);
    d = getStimvol(v,varname,'taskNum',taskNum,'phaseNum',phaseNum,...
        'segmentNum',segmentNum);
    assert(~isempty(d),'no events found');
    d_ind = 1:numel(d);
    if ~isempty(ignorelabels)
        d(d_ind==ignorelabels) = [];
        d_ind(d_ind==ignorelabels) = [];
    end
    conind = arrayfun(@(x)ones(size(d{x}))*x,1:numel(d),'uniformoutput',0);
    % -1 because mrTools is 1 indexed. (time 0 is scan 1)
    X{s} = convolveonsets(cat(2,d{:})'-1,cat(2,conind{:})',...
        epivol.frameperiod,viewGet(v,'nframes'),'dur',dur,...
        'hrffun',hrffun);
end
names = mat2strcell(d_ind,[varname '_%02d']);
designvol = Volume(cat(1,X{:}),...
    'metasamples',struct('chunks',epivol.meta.samples.chunks),...
    'metafeatures',struct('labels',{names}),...
    'frameperiod',epivol.frameperiod);

if madeview
    deleteView(v);
end
