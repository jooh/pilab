% Convert the data in a SPM.mat (path or loaded mat file) to epivol and
% designvol instances for further pilab processing.
%
% Named, optional inputs
% mask
% ignorelabels
%
% [epivol,designvol] = spm2vol(SPM,[varargin])
function [epivol,designvol] = spm2vol(SPM,varargin)

getArgs(varargin,{'mask',[],'ignorelabels',[]});

if ischar(SPM)
    % assume path to mat file
    SPM = loadbetter(SPM);
end

% construct epivol
epivol = SPMVolume(SPM,mask);

% construct designvol
X = SPM.xX.X;
nrun = numel(SPM.Sess);
names = arrayfun(@(x)x.name{1},SPM.Sess(1).U,...
  'uniformoutput',false);
[names,nameind] = setdiff(names,ignorelabels,'stable');

chunks = [];
for r = 1:nrun
    % task regressors only
    colind = SPM.Sess(r).col(1:length(SPM.Sess(r).U));
    % non-nuisance only
    colind = colind(nameind);
    % check that regressor names match across runs
    runnames = arrayfun(@(x)x.name{1},SPM.Sess(r).U,'uniformoutput',0);
    assert(isequal(names,runnames(nameind)),['regressor names must be ' ...
        'identical across runs']);
    runx{r} = X(SPM.Sess(r).row,colind);
    chunks = [chunks; ones(SPM.nscan(r),1)*r];
end
assert(isequal(chunks,epivol.meta.samples.chunks),['mismatched chunks ' ...
    'between epivol and designvol']);

designvol = Volume(vertcat(runx{:}),'metasamples',...
    struct('chunks',chunks),'metafeatures',struct('labels',{names}),...
    'frameperiod',SPM.xY.RT);
