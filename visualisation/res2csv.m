% export an roidata-style struct (e.g. outputs from roidata_rfx) to a CSV
% file with reasonable Excel/SPSS compatibility. Wrapper for table2csv.
%
% INPUT         DEFAULT         DESCRIPTION
% res           -               struct returned by roidata_* function
% filename      -               full path to output file
%
% NAMED INPUT       DEFAULT         DESCRIPTION
% target            'mean'          res field to export
% roitarget         res.cols_roi    names of ROIs to output (default all)
% contrasttarget    res.rows_contrast as roitarget
% subjecttarget     res.z_subject   as roitarget
% dotranspose       false           optionally flip the dims
% precision         []              if present, round to precision decimal points
% precfun           @round          function for rounding (ceil is good for p)
%
% res2csv(res,filename,[varargin])
function res2csv(res,filename,varargin)

getArgs(varargin,{'target','mean',...
    'precision',[],'roitarget',res.cols_roi,'contrasttarget',...
    res.rows_contrast,'subjecttarget',[],'dotranspose',false,...
    'precfun',[]);

% figure out which part of the matrix we're plotting
[~,~,roiind] = intersect(roitarget,res.cols_roi,'stable');
[~,~,conind] = intersect(contrasttarget,res.rows_contrast,'stable');
nroi = numel(roiind);
ncon = numel(conind);
assert(nroi==numel(roitarget),'did not find some roitarget entries');
assert(ncon==numel(contrasttarget),...
    'did not find some contrasttarget entries');

isgroupres = ndims(res.(target)) > 2;
if isgroupres
    if isempty(subjecttarget)
        % can only parse the subject input here since this field is missing
        % in many cases
        subjecttarget = res.z_subject;
        [~,~,subind] = intersect(subjecttarget,res.z_subject,'stable');
        nsub = numel(subind);
        assert(nsub==numel(subjecttarget),...
            'did not find some subjecttarget entries');
    end
    % need to combine the roiind and conind to a combined column axis
    out.x = cell(1,nroi * ncon);
    out.y = subjecttarget;
    % there is probably a more elegant way to achieve this
    out.data = NaN([length(out.y),length(out.x)]);
    n = 0;
    for r = 1:numel(roiind)
        for c = 1:numel(conind)
            n = n+1;
            out.x{n} = [roitarget{r},'_',contrasttarget{c}];
            out.data(:,n) = res.(target)(conind(c),...
                roiind(r),:);
        end
    end
    col1 = 'subject';
else
    % easy peasy
    out.data = res.(target)(conind,roiind,1);
    out.y = contrasttarget;
    out.x = roitarget;
    col1 = 'contrast';
end

if dotranspose
    out.data = out.data';
    [out.x,out.y] = deal(out.y,out.x);
end

% off to master function
table2csv(out.data,filename,'precision',precision,'precfun',precfun,...
    'collabels',out.x,'rowlabels',out.y);
