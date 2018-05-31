% upcast a split data rdm (n/2 by n/2 with no 0 diagonal and nrdm in 3rd
% dim) to a full rdm (n by n by nrdm with a 0 diagonal and the scalar
% replacewith (default NaN) in upper left and bottom right quadrants).
%
% NB! we assume that you the split data RDM is the the upper right quadrant
% of the original RDM, ie rdm(1:n,n+1:end). If this assumption is incorrect
% the split RDM will be inserted transposed with regard to the original
% configuration.
%
% fullrdm = upcastsplitrdm(splitrdm,[replacewith]);
function fullrdm = upcastsplitrdm(splitrdm,replacewith);

if ieNotDefined('replacewith')
    replacewith = NaN;
end
[nhalf,nhalf,nrdm] = size(splitrdm);
nfull = nhalf*2;

if isnan(replacewith)
    fullrdm = NaN([nfull nfull nrdm],class(splitrdm));
else
    fullrdm = ones([nfull nfull nrdm],class(splitrdm))*replacewith;
end

fullrdm(repmat(logical(eye(nfull)),[1 1 nrdm])) = 0;
fullrdm(1:nhalf,nhalf+1:end,:) = splitrdm;
fullrdm(nhalf+1:end,1:nhalf,:) = permute(splitrdm,[2 1 3]);
