% return a contrasts struct for roidata_rfx containing all pairwise
% contrasts between the conditions in the cell array connames.
%
% contrasts = roidata_allpairwisecontrasts(connames)
function contrasts = roidata_allpairwisecontrasts(connames)

ncon = numel(connames);
conmat = allpairwisecontrasts(ncon);
ncontrast = size(conmat,1);

for c = 1:ncontrast
    thiscon = conmat(c,:);
    conp = connames(thiscon==1);
    conn = connames(thiscon==-1);
    assert(numel(conp)==1 && numel(conn)==1,'contrast indexing failed');
    contrasts(c).name = sprintf('contrast_%s_o_%s',conp{1},conn{1});
    contrasts(c).conplus = conp;
    contrasts(c).conminus = conn;
    contrasts(c).tail = 'both';
end
