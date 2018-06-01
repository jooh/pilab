% resample the input rows/coluns of the input rdm(s) with replacement. Mask
% diagonal values that get resampled to off-diagonal positions with NaNs.
%
% rdm = resamplerdm(rdm,ind)
function rdm = resamplerdm(rdm,ind)

rdm = asrdmmat(rdm);
[ncon,~,nroi] = size(rdm);

% try out resampling and note where diagonal entries end up getting shifted off
% center
referencerdm = zerodiagonal(true([ncon,ncon]));
referencerdm = referencerdm(ind,ind);
matmask = repmat(zerodiagonal(referencerdm==0),[1 1 nroi]);

% resample
rdm = rdm(ind,ind,:);
rdm(matmask) = NaN;
