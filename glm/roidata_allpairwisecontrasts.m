% return a contrasts struct for roidata_rfx containing all pairwise
% contrasts between the conditions in the cell array connames.
%
% contrasts = roidata_allpairwisecontrasts(connames,convecmode)
function contrasts = roidata_allpairwisecontrasts(connames,convecmode)

if ieNotDefined('convecmode')
    convecmode = false;
end

ncon = numel(connames);
conmat = allpairwisecontrasts(ncon);
ncontrast = size(conmat,1);

for c = 1:ncontrast
    thiscon = conmat(c,:);
    conp = connames(thiscon==1);
    conn = connames(thiscon==-1);
    assert(numel(conp)==1 && numel(conn)==1,'contrast indexing failed');
    contrasts(c).name = sprintf('contrast_%s_o_%s',conp{1},conn{1});
    if convecmode
        contrasts(c).convec = conmat(c,:);
    else
        contrasts(c).conplus = conp;
        contrasts(c).conminus = conn;
    end
    contrasts(c).tail = 'both';
end
