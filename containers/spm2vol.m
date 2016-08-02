% Convert the data in a SPM.mat (char if path or loaded mat file) to epivol and
% designvol instances for further pilab processing.
%
% OPTIONAL INPUTS
% mask: see MriVolume
% ignorelabels: cell array of conditions to ignore
% blockmode=false: block concatenate runs (false means vertical concat).
%
% [epivol,designvol] = spm2vol(SPM,[varargin])
function [epivol,designvol] = spm2vol(SPM,varargin)

getArgs(varargin,{'mask',[],'ignorelabels',[],'blockmode',false});

if ischar(SPM)
    % assume path to mat file
    SPM = loadbetter(SPM);
end

% construct epivol
epivol = SPMVolume(SPM,mask);

% construct designvol
X = SPM.xX.X;
nrun = numel(SPM.Sess);

for r = 1:nrun
    allcolind = SPM.Sess(r).col(1:length(SPM.Sess(r).U));
    % non-nuisance only
    allnames = arrayfun(@(x)x.name{1},SPM.Sess(r).U,...
      'uniformoutput',false);
    [names{r},nameind] = setdiff(allnames,ignorelabels,'stable');
    colind{r} = allcolind(nameind);
    rowind{r} = SPM.Sess(r).row;
    chunks{r} = ones(SPM.nscan(r),1) * r;
end
chunks = vertcat(chunks{:});

if blockmode
    % block concatenated design matrix
    X = X(horzcat(rowind{:}),horzcat(colind{:}));
    names = vertcat(names{:});
else
    % vertical concatenation of design matrix
    X = arrayfun(@(x)X(rowind{x},colind{x}),1:nrun,'uniformoutput',0);
    X = vertcat(X{:});
    % check that names match
    assert(isequal(names{1},names{:}),['regressor names must be ' ...
        'identical across runs']);
    names = names{1};
end
assert(isequal(chunks,epivol.meta.samples.chunks),['mismatched chunks ' ...
    'between epivol and designvol']);
designvol = Volume(X,'metasamples',struct('chunks',chunks),...
    'metafeatures',struct('labels',{names}),'frameperiod',SPM.xY.RT);
