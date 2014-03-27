% serial version of runrois_spmd. Mainly useful for debugging ROIProcessor
% instances. See runrois_spmd for documentation.
%
% result = runrois(roimat,designmat,epimat,chunks,processor,minn,nreturn)
function result = runrois(roimat,designmat,epimat,chunks,processor,minn,nreturn)
if ischar(processor)
    processor = loadbetter(processor);
end
if isstruct(processor)
    processor = struct2obj(processor);
end
validvox = full(roimat~=0);
nvalid = sum(validvox,2);
nroi = size(validvox,1);
goodres = nvalid>0 & nvalid >=minn;
result = cell(nroi,nreturn);
retind = 1:nreturn;
for b = 1:nroi
    if ~goodres(b)
        continue
    end
    [result{b,retind}] = call(processor,designmat,...
        epimat(:,validvox(b,:)),chunks);
end
result(~goodres,:) = {NaN(size(result{find(...
    goodres,1,'first')}))};
