% run a processor on each of a set of ROIs (rows in roimat) using spmd
% parallelisation. Tries to be clever about only passing data to workers
% once, and only passing the bits of datamat that each worker needs. With
% reasonably big datasets (e.g. fMRI size - a few GB) this leads to a
% speedup over parfor and less risk of running out of memory on the
% workers.
%
% INPUTS:
% roimat: nroi by nfeatures (can be sparse)
% datamat: nsamples by nfeatures
% processor: function/method/struct that get's fevaled (see below)
% nreturn: number of returns from processor
% minn: minimum ROI size to run (smaller get NaNed out)
% [varargin]: any varargin get passed to processor (see below)
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
% However the processor is constructed, it must have a 'call' method/field
% which supports the following syntax:
% [batchresult{b,retind}] = call(processor,...
%   Cbatchdata(:,Cbatchfeat(b,:)),varargin{:});
%
% Note that if you are uncomfortable with object classes you can probably
% make your processor a struct with a function handle in the 'call' field.
% Matlab has trouble saving structs with function handles too, but you can
% use obj2struct to convert such structs to a format that can be passed in
% spmd (see above).
%
function result = runrois_spmd(roimat,datamat,processor,nreturn,minn,varargin)

retind = 1:nreturn;

% (nb transpose because distributed operates along last dim
% (columns)
Dvalidfeat = distributed(roimat');
% find the mask for each batch and restrict each workers roimat to these
spmd
    Cbatchfeat = full(getLocalPart(Dvalidfeat)'~=0);
    Cbatchmask = any(Cbatchfeat,1);
    Cbatchfeat = Cbatchfeat(:,Cbatchmask);
    Cnvalid = sum(Cbatchfeat,2);
    Cnroi = size(Cbatchfeat,1);
    Cbatchgoodres = (Cnvalid~=0 & Cnvalid>=minn);
end
clear Dvalidfeat;

Cbatchdata = Composite;
% this is not spmd to avoid passing the full datamat to workers
for batch = 1:numel(Cbatchdata)
    % retrieve the (small) Cbatchmask from worker and use to
    % index datamat
    Cbatchdata{batch} = datamat(:,Cbatchmask{batch});
end
% don't need this anymore
clear datamat;
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
    % get those feat (and transpose back)
    batchresult = cell(Cnroi,nreturn);
    for b = 1:Cnroi
        % check for invalid ROIs here
        if ~Cbatchgoodres(b)
            % empty or too small roi
            continue
        end
        [batchresult{b,retind}] = call(Cbatchprocessor,...
            Cbatchdata(:,Cbatchfeat(b,:)),varargin{:});
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
