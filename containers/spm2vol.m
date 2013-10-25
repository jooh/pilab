% Convert the data in a SPM.mat (path or loaded mat file) to epivol and
% designvol instances for further pilab processing.
%
% Named, optional inputs (see spm2glmdenoise for descriptions)
% mask
% ignorelabels
% assumeconvolved
% includeextras
%
% [epivol,designvol] = spm2vol(SPM,[varargin])
function [epivol,designvol] = spm2vol(SPM,varargin)

getArgs(varargin,{'mask',[],'ignorelabels',[],'assumeconvolved',[],...
    'includeextras',[]});

if ischar(SPM)
    % assume path to mat file
    SPM = loadbetter(SPM);
end

header = spm_vol(SPM.xY.P(1,:));

% this is a bit quirky since we won't glmdenoise but this function
% does a lot of heavy lifting that we don't want to replicate here
[epi,models,dur,mask,labels] = spm2glmdenoise(SPM,mask,ignorelabels,...
    assumeconvolved,includeextras);
% glmdenoise stores epis with samples in columns and features in
% rows
epi = cellfun(@(x)x',epi,'uniformoutput',false);
% make one matrix
epi = vertcat(epi{:});

nchunks = numel(models);
% set up chunks
chunks = cell2mat(arrayfun(@(x)ones(size(models{x},1),1)*x,...
    (1:nchunks)','uniformoutput',false));

% construct volume instances
epivol = MriVolume(epi,mask,'header',header,...
    'metasamples',struct('chunks',chunks,'volnames',...
    {{SPM.xY.VY.fname}'}),'frameperiod',SPM.xY.RT);
designvol = BaseVolume(vertcat(models{:}),'metafeatures',struct(...
    'labels',{labels}),'metasamples',struct('chunks',chunks),...
    'frameperiod',SPM.xY.RT);
