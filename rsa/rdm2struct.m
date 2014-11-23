% convert an RDM in matrix form (stacked in third dimension) to struct
% array.
%
% rdmstruct = rdm2struct(rdm,names)
function rdmstruct = rdm2struct(rdm,names)

[nc,~,nrdm] = size(rdm);

rdc = mat2cell(rdm,nc,nc,ones(nrdm,1));

rdmstruct = struct('name',names(:)','RDM',rdc(:)');
