% maskepivol = mask2meantimecourse(rois,epivol)
function maskepivol = mask2meantimecourse(rois,epivol)

% make sure we have the same ROIs and voxels across splits
[rois,epivol] = intersectvols(rois,epivol);
% now that the ROIs and voxels are in register this should reduce
% memory use considerably
validvox = any(rois.data~=0,1);
% at the moment we assume this for detecting searchlight maps. So
% need to crash if this is not met
if rois.nsamples == rois.nfeatures
    assert(all(validvox),'unused voxels in searchlight map');
else
    rois = rois(:,validvox);
    epivol = epivol(:,validvox);
end
assert(isequal(epivol.xyz,rois.xyz),...
    'mismatched roi/epivol');

% average the second input (ie, the EPIs from an ROI) over the second dim
fp = FunProcessor(@(data)mean(data,2),[],1);
roip = ROIProcessor(rois,fp);
meanres = call(roip,epivol.data);
maskepivol = result2roivol(roip,meanres,false,epivol.meta.samples);
