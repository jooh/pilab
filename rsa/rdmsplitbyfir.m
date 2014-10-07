% [outrdm,outnames] = rdmsplitbyfir(rdm,splitrdm)
function [outrdm,outnames] = rdmsplitbyfir(rdm,splitrdm)

rdmat = asrdmmat(rdm);
splitmat = asrdmmat(splitrdm);
nrdm = size(rdmat,3);
nsplit = size(splitmat,3);

if isstruct(rdm)
    names = {rdm.name};
    wasstruct = true;
else
    names = mat2strcell(1:nrdm,'rdm%02d');
    wasstruct = false;
end

if isstruct(splitrdm)
    splitnames = {splitrdm.name};
else
    splitnames = mat2strcell(1:nsplit,'split%02d');
end

for n = 1:nsplit
    % indices to use
    inds = splitmat(:,:,n) == 0;
    thisrdm = rdmat;
    thisrdm(repmat(inds,[1 1 nrdm])) = 0;
    outrdm{n} = thisrdm;
    outnames{n} = cellfun(@(thisname)[splitnames{n} '_' thisname],...
        names,'uniformoutput',false);
end
outrdm = cat(3,outrdm{:});
outnames = cat(2,outnames{:});

badtest = std(asrdmvec(outrdm))==0;
if any(badtest)
    logstr(...
        'removed %d RDMs which were invariant after splitting\n',...
        sum(badtest));
    outrdm(:,:,badtest) = [];
    outnames(badtest) = [];
end

if wasstruct
    outrdm = rdm2struct(outrdm,outnames);
end
