% Write out the RDM as a CSV file. Any varargin are passed to table2csv.
%
% INPUT     DEFAULT     DESCRIPTION
% rdm       -           something where isrdm returns true
% outpath   -           write CSV to this path
% names     {}          used as both row and column labels
% varargin  {}          any additional inputs are passed to table2csv.
%                           zlabels is especially useful for writing out
%                           multiple RDMs into the same file.
%
% rdm2csv(rdm,outpath,[names],[varargin])
function rdm2csv(rdm,outpath,names,varargin)

if ieNotDefined('names')
    names = {};
end

% make sure mat from here on
rdm = asrdmmat(rdm);

% NB need to insert extra column label to account for rowlabels
table2csv(rdm,outpath,'rowlabels',ascol(names),...
    'collabels',[{''},asrow(names)],varargin{:});
