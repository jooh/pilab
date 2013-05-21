% generate leave-one-out estimates of RDM reproducibility (spearman rho)
% for each entered dissimilarity volume.
%
% INPUT:
% comma-separated list of dissimilarity vols, e.g. disvolcell{:}. See
% roidata2rdmvol.
%
% OUTPUT:
% rbyroi: mean rho for each ROI (fisher transformed before averaging, then
%   back to rho)
% rbysplit: mean rho for each split (fisher transformed as above)
% allr: nrdm by nroi matrix of rho
%
% [rbyroi,rbysplit,allr] = rdmreproducibility(varargin)
function [rbyroi,rbysplit,allr] = rdmreproducibility(varargin)

nrdm = nargin;
assert(nrdm>1,'must enter at least 2 rdms');

ndis = varargin{1}.nsamples;
nroi = varargin{1}.nfeatures;

% verify that the disvols are homogeneous
names = cellfun(@(disvol)disvol.meta.features.names,varargin,...
    'uniformoutput',false);
assert(isequal(names{:}),'mismatched ROI names across disvols')
ndisall = cellfun(@(disvol)disvol.nsamples,varargin,'uniformoutput',true);
assert(all(ndisall==ndis),'mismatched RDMs across disvols')
nroiall = cellfun(@(disvol)disvol.nfeatures,varargin,'uniformoutput',true);
assert(all(nroiall==nroi),'mismatched ROIs across disvols')

% iterate over ROIs
allr = NaN([nrdm nroi]);
for roi = 1:nroi
    % pull RDMs into matrix
    rdms = cellfun(@(disvol)disvol.data(:,roi),varargin,'uniformoutput',...
        false);
    rdmat = cell2mat(rdms);
    % generate leave-one-out estimates for each fold
    for sp = 1:nrdm
        train = mean(rdmat(:,(1:nrdm)~=sp),2);
        allr(sp,roi) = corr(train,rdmat(:,sp),'type','spearman');
    end
end
% make sure any perfect rs don't produce inf results by rounding down to 4
% decimal places (Z of around 5)
toostrong = abs(allr)==1;
allr(toostrong) = allr(toostrong) * .9999;
allfish = atanh(allr);
rbyroi = tanh(mean(allfish,1));
rbysplit = tanh(mean(allfish,2));
