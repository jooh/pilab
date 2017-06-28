% export experimental events from SPM.mat file to BIDS-compliant TSV file.
%
% INPUT         DEFAULT     DESCRIPTION
% SPM           -           SPM struct or path to one
% outroot       -           leading path, e.g. '/path/to/sub-01_task-main'. We
%                               append 'run-x_events.tsv' to make the final
%                               output file name per run.
% ignorelabels  {}          regressors to skip
% selectruns    1:nrun      runs to include (note that we renumber, so
%                               selectruns=[1,3,5] will get written out as
%                               run_01, run_02, run_03)
%
% 2017-05-23 J Carlin
%
% spmevents2bids(SPM,outroot,[ignorelabels],[selectruns])
function spmevents2bids(SPM,outroot,ignorelabels,selectruns)

if ischar(SPM)
    SPM = load(SPM);
    SPM = SPM.SPM;
end

if ~exist('ignorelabels','var') || isempty(ignorelabels)
    ignorelabels = {};
end
if ~iscell(ignorelabels)
    ignorelabels = {ignorelabels};
end

nrun = numel(SPM.Sess);
if ~exist('selectruns','var') || isempty(selectruns)
    selectruns = 1:nrun;
end

switch lower(SPM.xBF.UNITS)
    case 'scans'
        unitconv = @(x)x ./ SPM.xBF.RT;
    case 'secs'
        unitconv = @(x)x;
    otherwise
        error('unknown SPM.xBF.units: %s',SPM.xBF.units);
end

rc = 0;
for r = selectruns(:)'
    rc = rc+1;
    onsets = [];
    names = {};
    dur = [];
    for thisu = SPM.Sess(r).U(:)'
        if any(ismember(ignorelabels,thisu.name{1}))
            continue;
        end
        names = [names repmat(thisu.name,[1 numel(thisu.ons)])];
        onsets = [onsets unitconv(thisu.ons)];
        dur = [dur unitconv(thisu.dur)];
    end
    % sort by time
    [onsets,sortind] = sort(onsets);
    names = names(sortind);
    dur = dur(sortind);
    outfile = sprintf('%s_run-%02d_events.tsv',outroot,rc);
    writebidsevents(outfile,onsets,dur,names);
end
