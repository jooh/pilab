% generate nboot samples with replacement of the indices in 1:ncon. 
%
% To avoid problems with invariant features, we drop any bootstraps where
% the same single sample is drawn n times.
%
% bootinds = bootindices(ncon,nboot)
function bootinds = bootindices(ncon,nboot)

npossible = (ncon^ncon) - ncon;
if nboot > npossible
    warning(...
    'requested %d bootstraps but only %d are possible',nboot,npossible);
    nboot = npossible;
end

if ncon<=8
    % just compute all the bootstraps and sample from this (memory
    % expensive but fast)
    bootinds = npermutek(1:ncon,ncon);
    % drop the non-unique draws
    nonunique = all(bootinds==repmat(bootinds(:,1),[1 ncon]),2);
    bootinds(nonunique,:) = [];
    bootinds = bootinds(randperm(npossible),:);
    bootinds = bootinds(1:nboot,:);
else
    % draw random, non-repeating boots 
    done = false;
    niter = 0;
    while ~done
        bootinds = nuniquereturns(@()ceil(rand(1,ncon)*ncon),nboot,...
            [1 ncon]);
        bootinds = vertcat(bootinds{:});
        % it's exceedingly unlikely in most scenarios but just to be safe
        % and prevent weird, unpredictable crashes later on
        nonunique = all(bootinds==repmat(bootinds(:,1),[1 ncon]),2);
        if ~any(nonunique)
            done = true;
        end
        niter = niter + 1;
        assert(niter<1e3,'iteration limit exceeded');
    end
end
