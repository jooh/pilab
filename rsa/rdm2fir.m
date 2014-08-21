% [firdm,outnames] = rdm2fir(rdm)
function [firdm,outnames] = rdm2fir(rdm)

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

for f = 1:nfir
    newvec = zeros([ndis 1],class(rdvec));
    newvec(rdvec==u(f)) = 1;
    firdm(f) = struct('name',sprintf('%s_f%02.0f',name,u(f)),...
        'RDM',vec2rdm(newvec));
end

outnames = {firdm.name};
if ~wasstruct
    firdm = asrdmmat(firdm);
end
