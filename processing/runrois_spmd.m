% run a processor on some set of ROIs (rows in roimat) using spmd
% parallelisation. Tries to be clever about only passing data to workers
% once, and only passing the bits of epimat that each worker needs.
%
% INPUTS:
% roimat: nroi by nfeatures (can be sparse)
% designmat: nsamples by nconditions
% epimat: nsamples by nfeatures
% chunks: nsamples by 1 vector
% processor: function/method/struct that get's fevaled (see below)
% minn: minimum ROI size to run (smaller get NaNed out)
% nreturn: number of returns from processor
%
% OUTPUTS:
% result: nsamples by nreturn cell array
%
% PROCESSOR HANDLING
% Matlab has serious problems with saving handle objects (and spmd
% basically uses save internally to pass variables).  We therefore support
% a few hacky ways to pass the processor: 
% a) the processor can be a char, which gets loaded on each worker 
% b) the processor (or the result of loading the above) can be a struct,
%   which gets converted to object instance by struct2obj
%
% However the processor is constructed, it then gets called like so:
% [batchresult{b,retind}] = call(processor,designmat,...
%   Cbatchepi(:,Cbatchvox(b,:)),chunks);
%
% Note that if you are uncomfortable with object classes you can probably
% make your processor a struct with a function handle in the 'call' field.
% Matlab has trouble saving structs with function handles too, but you can
% use obj2struct to convert such structs to a format that can be passed in
% spmd (see above).
%
% result = runrois_spmd(roimat,designmat,epimat,chunks,processor,...
        % minn,nreturn)
function result = runrois_spmd(roimat,designmat,epimat,chunks,processor,...
        minn,nreturn)
retind = 1:nreturn;

% (nb transpose because distributed operates along last dim
% (columns)
Dvalidvox = distributed(roimat');
% find the mask for each batch and restrict each workers roimat to these
spmd
    Cbatchvox = full(getLocalPart(Dvalidvox)'~=0);
    Cbatchmask = any(Cbatchvox,1);
    Cbatchvox = Cbatchvox(:,Cbatchmask);
    Cnvalid = sum(Cbatchvox,2);
    Cnroi = size(Cbatchvox,1);
    Cbatchgoodres = (Cnvalid~=0 & Cnvalid>=minn);
end
clear Dvalidvox;

Cbatchepi = Composite;
% this is not spmd to avoid passing the full epimat to workers
for batch = 1:numel(Cbatchepi)
    % retrieve the (small) Cbatchmask from worker and use to
    % index epimat
    Cbatchepi{batch} = epimat(:,Cbatchmask{batch});
end
% don't need this anymore
clear epimat;
clear Cbatchmask;

spmd
    if ischar(processor)
        Cbatchprocessor = loadbetter(processor);
        if isstruct(Cbatchprocessor)
            Cbatchprocessor = struct2obj(Cbatchprocessor);
        end
    elseif isstruct(processor)
        Cbatchprocessor = struct2obj(processor);
    else
        Cbatchprocessor = processor;
    end
    % get those voxels (and transpose back)
    batchresult = cell(Cnroi,nreturn);
    for b = 1:Cnroi
        % check for invalid ROIs here
        if ~Cbatchgoodres(b)
            % empty or too small roi
            continue
        end
        [batchresult{b,retind}] = call(Cbatchprocessor,designmat,...
            Cbatchepi(:,Cbatchvox(b,:)),chunks);
    end
end

result = cat(1,batchresult{:});
goodres = cat(1,Cbatchgoodres{:});
clear batchresult
clear C*
% insure that we plug in appropriate nans
result(~goodres,:) = {NaN(size(result{find(...
    goodres,1,'first')}))};
% there are a few more variables that could be cleared but I
% think this will do for now
