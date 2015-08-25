% Write out the RDM as a CSV file. Any varargin are passed to table2csv.
%
% rdm2csv(rdm,outpath,[names],[varargin])
function rdm2csv(rdm,outpath,names,varargin)

if ieNotDefined('names')
    names = [];
end

if isstruct(rdm) & isempty(names)
    names = {rdm.name};
end

% make sure mat from here on
rdm = asrdmmat(rdm);

[nc,~,nrdm] = size(rdm);
assert(nrdm==1,'only 1 input RDM can be entered per call');

% NB need to insert extra column label to account for rowlabels
table2csv(rdm,outpath,'rowlabels',names,'collabels',[{''},names(:)'],...
    varargin{:});
