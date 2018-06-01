% split data RDM from a 2 run master RDM. TODO: support for combining
% arbitrary n splits to one RDM.
% rdm = splitrdm(rdm)
function rdm = splitrdm(rdm)

n = size(rdm,1);
midp = n/2;
rows = (midp+1):n;
cols = 1:midp;
rdm = rdm(rows,cols);
