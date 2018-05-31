% serial version of runrois_spmd (see this function for documentation).
%
% This code should be functionally equivalent to runrois_spmd but will a)
% be easier to use when debugging processors, b) run a lot faster when you
% already have parallelisation in your code -  (since Matlab only
% parallelises the outermost level and all those unused parallel toolbox
% bits in runrois_spmd still come with overhead).
%
% result = runrois_serial(roimat,datamat,processor,nreturn,minn,[varargin])
function result = runrois_serial(roimat,datamat,processor,nreturn,minn,varargin)

if ischar(processor)
    processor = loadbetter(processor);
end
if isstruct(processor)
    processor = struct2obj(processor);
end

validfeat = full(roimat~=0);
nvalid = sum(validfeat,2);
nroi = size(validfeat,1);
goodres = nvalid>0 & nvalid >=minn;
result = cell(nroi,nreturn);
retind = 1:nreturn;
for b = 1:nroi
    if ~goodres(b)
        continue
    end
    [result{b,retind}] = call(processor,...
        datamat(:,validfeat(b,:)),varargin{:});
end
result(~goodres,:) = {NaN(size(result{find(...
    goodres,1,'first')}))};
