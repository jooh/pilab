% split the input rdm to one output firdm per unique value in the input.
% Numformat determines how we name the output rdms. By default output RDMs
% are scored using 1 where the unique value was in the input and 0
% elsewhere. If keepvalues is true we instead plug in the actual
% dissimilarity.
% [firdm,outnames] = rdm2fir(rdm,[numformat='%02.0f'],[keepvalues=0])
function [firdm,outnames] = rdm2fir(rdm,numformat,keepvalues)

if ieNotDefined('numformat')
    numformat = '%02.0f';
end

if ieNotDefined('keepvalues')
    keepvalues = 0;
end

name = 'firdm';
wasstruct = false;
if isstruct(rdm)
    name = rdm(1).name;
    wasstruct = true;
end

rdvec = asrdmvec(rdm);
[ndis,nrdm] = size(rdvec);
assert(nrdm==1,'only one input RDM is supported');
assert(all(std(rdvec)>0),'invariant rdms detected');

u = unique(rdvec);
u(isnan(u)) = [];
nfir = numel(u);

outvals = ones(nfir,1);
if keepvalues == 1
    % use the actual dissimilarities instead.
    outvals = u;
end

for f = 1:nfir
    newvec = zeros([ndis 1],class(rdvec));
    newvec(rdvec==u(f)) = outvals(f);
    firdm(f) = struct('name',sprintf(['%s_' numformat],name,u(f)),...
        'RDM',vec2rdm(newvec));
end

outnames = {firdm.name};
if ~wasstruct
    firdm = asrdmmat(firdm);
end
