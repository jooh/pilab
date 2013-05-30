% given a set of ROIs (searchlight spheres or masks) and a set of response
% pattern estimates, make a volume of (vectorised) dissimilarities
% according to distancemetric (default corr).
% disvol = roidata2rdmvol(rois,data,[distancemetric])
function disvol = roidata2rdmvol(rois,vol,distancemetric)

if ieNotDefined('distancemetric')
    distancemetric = 'corr';
end

% preallocate result (dissimilarity by roi)
dissimilarities = NaN([nchoosek(vol.desc.samples.nunique.labels,2) ...
    rois.nsamples]);

% now mask out any NaN features
nanmask = ~any(isnan(vol.data),1);
% if you haven't NaNed out all conditions for a given feature something is
% likely broken
assert(all(all(isnan(vol.data(:,~nanmask)))),'inconsistent NaN mask detected');
if ~all(nanmask)
    fprintf('removed %d NaN features from analysis.\n',sum(~nanmask));
end

% compute result
parfor n = 1:rois.nsamples
    % skip empty rois (these come out as NaN)
    if ~any(rois.data(n,:)&nanmask)
        continue
    end
    dissimilarities(:,n) = pdist(vol.data(:,full(rois.data(n,:)~=0)&nanmask),distancemetric);
end

% convert to volume - here it is a problem that the result may have
% different nfeatures than the mask (e.g. for ROI analysis or when we do
% not run all possible searchlight spheres)
if rois.nsamples == rois.nfeatures
    % simple case - assume that samples and features are in register
    disvol = MriVolume(dissimilarities,rois,'metafeatures',struct(...
        'names',{rois.meta.samples.names}));
else
    % complicated case - need to forget the mask and write out a mask-less
    % volume. But save coordinates of ROIs to enable sanity checks later
    coords = cell(1,rois.nsamples);
    nvox = NaN([1 rois.nsamples]);
    for c = 1:rois.nsamples
        % compute centre of mass for this ROI
        coords{c} = round(mean(rois.linind2coord(rois.linind(...
            rois.data(c,:)~=0)),2));
        nvox(c) = sum(rois.data(c,:)~=0);
    end
    % make a mask-less volume 
    disvol = MriVolume(dissimilarities,[],'metafeatures',struct(...
        'names',{rois.meta.samples.names'},'centreofmass',{coords},...
        'nfeatures',nvox),'header',rois.header);
end
