% export a roidata result (struct output from roidata_rsa, roidata_rfx etc) to
% CSV table.
%
% INPUT         DEFAULT         DESCRIPTION
% res           -               struct from roidata_rfx or similar
% outpath       -               path to for output CSV
% conind        1:ncon          numeric or logical indices into the rows of res.
%                                   We also support cell array inputs, in which
%                                   case we use intersect to find indices.
%                                   Finally, if you enter the char
%                                   'nocontrasts' we output everything that
%                                   doesn't begin with 'contrast_'.
% roiind        1:nroi          numeric or logical indices into the columns of
%                                   res.
% precision     []              passed on to table2csv - see that function
% varargin      'mean','sterr','ppara','n'  comma separated list of fields to
%                                               output nested inside each
%                                               condition.
%
% roidata2csv(res,outpath,[conind],[roiind],[precision],[varargin])
function roidata2csv(res,outpath,conind,roiind,precision,varargin)

% input parsing
targetfields = varargin;
if isempty(targetfields)
    targetfields = {'mean','sterr','ppara','n'};
end
ntarget = numel(targetfields);
[ncon,nroi,nsub] = size(res.(targetfields{1}));
assert(nsub==1,'no support for groupres at present')

if ieNotDefined('roiind')
    roiind = 1:nroi;
elseif iscell(roiind)
    % do intersection to find numerical index
    [~,~,roiind] = intersect(roiind,res.cols_roi,'stable');
end
if islogical(roiind)
    roiind = find(roiind);
end
assert(isnumeric(roiind),'roiind must be index')
nroi = numel(roiind);

if ieNotDefined('conind')
    conind = 1:ncon;
elseif iscell(conind)
    % do interssection to find numerical index
    [~,~,conind] = intersect(conind,res.rows_contrast,'stable');
elseif isstr(conind) && strcmp(lower(conind),'nocontrasts')
    % remove contrasts
    hits = strfindcell(res.rows_contrast,'contrast_');
    conind = setdiff(1:ncon,hits);
end
if islogical(conind)
    conind = find(conind);
end
assert(isnumeric(conind),'conind must be index')
ncon = numel(conind);

if ieNotDefined('precision')
    precision = [];
end

% extract data to plot. there are more efficient solutions here but hey ho
nrow = ncon * ntarget;
datamat = NaN([nrow,nroi]);
rowlabels = cell(nrow,2);
[rowlabels{:}] = deal('');
rc = 1;
for con = 1:ncon
    rowlabels{rc,1} = res.rows_contrast{conind(con)};
    for t = 1:ntarget
        rowlabels{rc,2} = targetfields{t};
        % get data
        datamat(rc,:) = res.(targetfields{t})(conind(con),roiind);
        rc = rc + 1;
    end
end

% and now perhaps we can just
table2csv(datamat,outpath,'precision',precision,'rowlabels',rowlabels,...
    'collabels',[{'',''} res.cols_roi(roiind)]);
