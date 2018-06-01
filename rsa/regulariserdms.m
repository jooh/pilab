% regularise (ie dimensionality reduce) datardms according to some model
% dissimilarities.  Finds all unique pairings of the regulariser's
% dissimilarities and replaces those entries in the datardms with the mean.
% This usually helps bring out structure in the RDM, but complicates
% hypothesis testing. As with any averaging, most tests will come out very
% similarly in terms of statistical significance with and without
% regularisation. So the main utility of this function is to aid
% visualisation of noisy RDMs.
%
% INPUTS:
% datardms: some RDM type (struct, vector, matrix)
% regulariser: some RDM type
%
% OUTPUTS:
% regrdms: vectorised RDMs
% udis: vector of averaged dissimilarities (non-RDM)
% ureg: vector of unique predictions in regularised (non-RDM)
%
% [regrdms,udis,ureg] = regulariserdms(datardms,regulariser)
function [regrdms,udis,ureg] = regulariserdms(rdms,regulariser)

% find the unique combinations of regulariser dissimilarities
regulariser = asrdmvec(regulariser);
% remove any NaN dissimilarities across all models.
ureg = uniquen(regulariser(any(~isnan(regulariser),2),:),'rows');
assert(~isempty(ureg),'no non-nan dissimilarities in regulariser');
nu = length(ureg);

rdms = asrdmvec(rdms);
[ndis,nrdm] = size(rdms);

regrdms = NaN([ndis,nrdm],class(rdms));
udis = NaN([nu nrdm]);

for u = 1:nu
    % find the dissimilarities that belong to this unique combination of
    % dissimilarities
    % (this is a bit ugly but checking for equality with NaNs is tricky)
    disind = arrayfun(@(x)isequaln(regulariser(x,:),ureg(u,:)),...
        [1:size(regulariser,1)]');
    for r = 1:nrdm
        udis(u,r) = mean(rdms(disind,r));
        regrdms(disind,r) = udis(u,r);
    end
end
