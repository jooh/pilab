% regularise (ie dimensionality reduce) datardms according to some model
% dissimilarities.  Finds all unique pairings of the regulariser's
% dissimilarities and replaces those entries in the datardms with the mean.
% This produces tons of structure in the RDM, but complicates hypothesis
% testing. 
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
ureg = unique(regulariser(any(~isnan(regulariser),2),:),'rows');
assert(~isempty(ureg),'no non-nan dissimilarities in regulariser');
nu = length(ureg);

rdms = asrdmvec(rdms);
[ndis,nrdm] = size(rdms);

regrdms = NaN([ndis,nrdm],class(rdms));
udis = NaN([nu nrdm]);

for u = 1:nu
    % find the dissimilarities that belong to this unique combination of
    % dissimilarities
    disind = all(bsxfun(@eq,regulariser,ureg(u,:)),2);
    for r = 1:nrdm
        udis(u,r) = mean(rdms(disind,r));
        regrdms(disind,r) = udis(u,r);
    end
end
