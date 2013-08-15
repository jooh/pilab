% generate nboot samples with replacement of the indices in 1:ncon. 
%
% bootinds = bootindices(ncon,nboot)
function bootinds = bootindices(ncon,nboot)

npossible = ncon^ncon;
if nboot > npossible
    warning(...
    'requested %d bootstraps but only %d are possible',nboot,npossible);
    nboot = npossible;
end

if ncon<=8
    % just compute all the bootstraps and sample from this (memory
    % expensive but fast)
    bootinds = npermutek(1:ncon,ncon);
    bootinds = bootinds(randperm(npossible),:);
    bootinds = bootinds(1:nboot,:);
else
    % draw random, non-repeating boots 
    bootinds = nuniquereturns(@()ceil(rand(1,ncon)*ncon),nboot,[1 ncon]);
    bootinds = vertcat(bootinds{:});
end
