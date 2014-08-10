% convert onsets in SPM struct to format supported by convolveonsets.
%
% INPUTS:
% SPM: loaded SPM.mat struct
% chunk: index in SPM.Sess to target
% [resortind]: custom ordering rather than follow regressor order in SPM
%
% OUTPUTS:
% onsets: vector of onsets in volume units
% conind: condition index for each onset (columns in SPM design)
% dur: duration of each onset in volume units
%
% [onsets,conind,dur] = spm2onsets(SPM,chunk,resortind)
function [onsets,conind,dur] = spm2onsets(SPM,chunk,resortind)

onsets = [];
conind = [];
dur = [];

ncon = length(SPM.Sess(chunk).U);
if ieNotDefined('resortind')
    resortind = 1:ncon;
end
assert(numel(resortind)==ncon,'resortind must be ncon length');

for con = 1:ncon
    thisons = SPM.Sess(chunk).U(con).ons;
    n = numel(thisons);
    thisdur = SPM.Sess(chunk).U(con).dur;
    if numel(thisdur)~=n
        if isscalar(thisdur)
            thisdur = repmat(thisdur,[n 1]);
        else
            error('duration does not match onsets');
        end
    end
    if strcmp(SPM.xBF.UNITS,'secs')
        thisons = thisons / SPM.xY.RT;
        thisdur = thisdur / SPM.xY.RT;
    end
    thisind = repmat(resortind(con),[n 1]);
    onsets = [onsets; thisons(:)];
    conind = [conind; thisind(:)];
    dur = [dur; thisdur(:)];
end

% sort by time
[onsets,sortind] = sort(onsets);
conind = conind(sortind);
dur = dur(sortind);
