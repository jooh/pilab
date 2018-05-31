% export experimental events to a BIDS-compliant TSV file.
%
% INPUT
% outfile
% onsets
% durations
% names
% extracols
%
% writebidsevents(outfile,onsets,durations,names,[extracols])
function writebidsevents(outfile,onsets,durations,names,extracols)

if ~exist('extracols','var') || isempty(extracols)
    extracols = emptystruct('name','value','format');
end

datacols = struct('name',{'onset','duration','trial_type'},...
    'value',{onsets,durations,names},'format',{'%.2f','%.2f','%s'});

datacols = [datacols extracols];

columns2file(outfile,datacols,'\t');
